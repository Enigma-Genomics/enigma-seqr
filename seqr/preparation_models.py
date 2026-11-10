"""Internal preparation idempotency state. Not exposed through clinical JSON APIs."""
from django.db import models


class PreparationLedger(models.Model):
    subject_key = models.CharField(max_length=128, unique=True)
    review_sequence = models.PositiveBigIntegerField(default=0)
    revision = models.CharField(max_length=64, blank=True, default='')
    request_id = models.CharField(max_length=64, blank=True, default='')
    project_guid = models.CharField(max_length=30, blank=True, default='')
    metadata_state_token = models.CharField(max_length=64, blank=True, default='')
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        app_label = 'seqr'


class PreparationOutbox(models.Model):
    ledger = models.ForeignKey(PreparationLedger, on_delete=models.PROTECT)
    revision = models.CharField(max_length=64)
    review_sequence = models.PositiveBigIntegerField()
    payload = models.JSONField()
    delivered_at = models.DateTimeField(null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        app_label = 'seqr'
        constraints = [
            models.UniqueConstraint(fields=['ledger', 'revision'], name='prep_outbox_revision_unique'),
            models.UniqueConstraint(fields=['ledger', 'review_sequence'], name='prep_outbox_sequence_unique'),
        ]
