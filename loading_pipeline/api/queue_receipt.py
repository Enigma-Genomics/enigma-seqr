"""Publish a validated request once on the runner's persistent queue volume."""
import fcntl
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path


class QueueFull(Exception):
    pass


class RequestConflict(ValueError):
    pass


def _sync(directory):
    fd = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _write(path, data):
    with path.open('x') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())


def publish_once(queue, key, payload, allocate_id, queue_path, capacity):
    if not re.fullmatch(r'[a-f0-9]{64}', key):
        raise ValueError('Invalid load request key')
    queue = Path(queue)
    receipts = queue.parent / 'loading_request_receipts'
    receipts.mkdir(exist_ok=True)
    canonical = json.dumps(json.loads(payload), sort_keys=True, separators=(',', ':'))
    fingerprint = hashlib.sha256(canonical.encode()).hexdigest()
    target = receipts / key
    with (receipts / '.publish.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not target.exists():
            if len(os.listdir(queue)) >= capacity:
                raise QueueFull
            run_id = allocate_id()
            with tempfile.TemporaryDirectory(dir=receipts) as scratch:
                staged = Path(scratch) / 'request'
                staged.mkdir()
                _write(staged / 'receipt.json', json.dumps({'run_id': run_id, 'fingerprint': fingerprint}))
                _write(staged / 'pending', canonical)
                _sync(staged)
                staged.rename(target)
                _sync(receipts)
        receipt = json.loads((target / 'receipt.json').read_text())
        if receipt['fingerprint'] != fingerprint:
            raise RequestConflict('This request key already identifies a different load')
        pending = target / 'pending'
        if pending.exists():
            # A single rename both publishes to Luigi and consumes the right to
            # publish. A replay after Luigi removes its file still sees the receipt.
            pending.rename(queue_path(receipt['run_id']))
            _sync(queue)
            _sync(target)
        return receipt['run_id']
