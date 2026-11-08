from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):
    dependencies = [('seqr', '0090_analysisgroup_workspace_name_and_more')]
    operations = [
        migrations.CreateModel(name='PreparationLedger', fields=[
            ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ('subject_key', models.CharField(max_length=128, unique=True)),
            ('review_sequence', models.PositiveBigIntegerField(default=0)),
            ('revision', models.CharField(blank=True, default='', max_length=64)),
            ('request_id', models.CharField(blank=True, default='', max_length=64)),
            ('project_guid', models.CharField(blank=True, default='', max_length=30)),
            ('metadata_state_token', models.CharField(blank=True, default='', max_length=64)),
            ('updated_at', models.DateTimeField(auto_now=True)),
        ]),
        migrations.CreateModel(name='PreparationOutbox', fields=[
            ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ('revision', models.CharField(max_length=64)),
            ('review_sequence', models.PositiveBigIntegerField()),
            ('payload', models.JSONField()),
            ('delivered_at', models.DateTimeField(blank=True, null=True)),
            ('created_at', models.DateTimeField(auto_now_add=True)),
            ('ledger', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='seqr.preparationledger')),
        ]),
        migrations.AddConstraint(model_name='preparationoutbox', constraint=models.UniqueConstraint(
            fields=('ledger', 'revision'), name='prep_outbox_revision_unique')),
        migrations.AddConstraint(model_name='preparationoutbox', constraint=models.UniqueConstraint(
            fields=('ledger', 'review_sequence'), name='prep_outbox_sequence_unique')),
    ]
