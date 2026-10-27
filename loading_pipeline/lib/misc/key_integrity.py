# Table expressions come from the trusted TableNameBuilder, never request input.
# ruff: noqa: S608

from collections.abc import Callable


class VariantKeyIntegrityError(ValueError):
    """Variant identity cannot be preserved by importing this run."""


def _reject_conflict(execute: Callable, query: str, reason: str) -> None:
    if execute(query):
        # Keep clinical variant identifiers out of exception and notification text.
        msg = (
            f'Variant key integrity check failed: {reason}. '
            'Import stopped before writes. Regenerate stale run outputs against a '
            'current variant snapshot, or repair the existing mapping first.'
        )
        raise VariantKeyIntegrityError(msg)


def validate_new_variant_keys(
    execute: Callable,
    source_table: str,
    lookup_table: str,
) -> None:
    """Reject stale key allocations before importing any variant tables.

    These checks do not reserve keys: callers must serialize writers from this
    check through the import, or use an atomic key-range allocator.
    """
    checks = [
        (
            f"""
            SELECT 1 FROM (
                SELECT 1 FROM {source_table}
                WHERE key IS NULL OR key < 0 OR key > 4294967295
                   OR variantId IS NULL OR variantId = ''
            ) LIMIT 1
            """,
            'invalid source key or variant identifier',
        ),
        (
            f"""
            SELECT 1 FROM (
                SELECT 1 FROM {source_table} GROUP BY key HAVING count() > 1
                UNION ALL
                SELECT 1 FROM {source_table} GROUP BY variantId HAVING count() > 1
            ) LIMIT 1
            """,
            'source mapping is not one-to-one',
        ),
        (
            f"""
            SELECT 1 FROM (
                SELECT 1 FROM {source_table} src
                INNER JOIN (
                    SELECT key, variantId FROM {lookup_table}
                    WHERE key IN (SELECT toUInt32(key) FROM {source_table})
                ) dst ON assumeNotNull(toUInt32(src.key)) = dst.key
                WHERE assumeNotNull(src.variantId) != dst.variantId
            ) LIMIT 1
            """,
            'source key is already assigned to a different variant',
        ),
        (
            f"""
            SELECT 1 FROM (
                SELECT 1 FROM {source_table} src
                INNER JOIN (
                    SELECT key, variantId FROM {lookup_table}
                    WHERE variantId IN (SELECT variantId FROM {source_table})
                ) dst ON assumeNotNull(src.variantId) = dst.variantId
                WHERE assumeNotNull(toUInt32(src.key)) != dst.key
            ) LIMIT 1
            """,
            'source variant is already assigned to a different key',
        ),
    ]
    for query, reason in checks:
        _reject_conflict(execute, query, reason)


def validate_entry_variant_keys(
    execute: Callable,
    source_table: str,
    lookup_table: str,
) -> None:
    """Reject entries whose identities are missing or ambiguous in the lookup."""
    checks = [
        (
            f"""
            SELECT 1 FROM (
                SELECT 1 FROM {source_table}
                WHERE variantId IS NULL OR variantId = ''
            ) LIMIT 1
            """,
            'entry has no variant identifier',
        ),
        (
            f"""
            WITH source_ids AS (SELECT DISTINCT variantId FROM {source_table})
            SELECT 1 FROM (
                SELECT 1 FROM source_ids src
                LEFT JOIN (
                    SELECT variantId, uniqExact(key) AS key_count
                    FROM {lookup_table}
                    WHERE variantId IN (SELECT variantId FROM source_ids)
                    GROUP BY variantId
                ) dst ON assumeNotNull(src.variantId) = dst.variantId
                WHERE ifNull(dst.key_count, 0) != 1
            ) LIMIT 1
            """,
            'entry variant does not resolve to exactly one key',
        ),
        (
            f"""
            SELECT 1 FROM (
                SELECT 1 FROM {lookup_table}
                WHERE key IN (
                    SELECT key FROM {lookup_table}
                    WHERE variantId IN (SELECT DISTINCT variantId FROM {source_table})
                )
                GROUP BY key HAVING uniqExact(variantId) > 1
            ) LIMIT 1
            """,
            'entry key is assigned to multiple variants',
        ),
    ]
    for query, reason in checks:
        _reject_conflict(execute, query, reason)
