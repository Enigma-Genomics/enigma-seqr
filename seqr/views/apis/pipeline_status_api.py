import json
import os
from datetime import datetime, timezone
from zoneinfo import ZoneInfo

from django.http import JsonResponse
from google.api_core.exceptions import GoogleAPICallError
from google.cloud import storage

from seqr.views.utils.permissions_utils import data_manager_required


@data_manager_required
def daily_sample_status(request):
    location = os.environ.get('PIPELINE_STATUS_SNAPSHOT_PREFIX')
    if not location:
        return JsonResponse({'enabled': False})
    now = datetime.now(timezone.utc)
    day = now.astimezone(ZoneInfo(os.environ.get('PIPELINE_STATUS_TIMEZONE', 'UTC'))).date().isoformat()
    bucket, prefix = location.removeprefix('gs://').rstrip('/').split('/', 1)
    try:
        blob = storage.Client().bucket(bucket).blob(f'{prefix}/{day}.json')
        snapshot = json.loads(blob.download_as_bytes(timeout=10))
        age = now - datetime.fromisoformat(snapshot['updated'])
        if snapshot['day'] != day:
            raise ValueError('Unexpected day')
        result = JsonResponse({**snapshot, 'enabled': True, 'stale': age.total_seconds() > 300})
    except (GoogleAPICallError, KeyError, ValueError):
        result = JsonResponse({'error': 'Daily sample status is temporarily unavailable.'}, status=503)
    result['Cache-Control'] = 'private, no-store'
    return result
