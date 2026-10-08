"""A radius on each campaign site, and the CampaignEvent projection.

Schema only. CampaignEvent is a cache of what finished pandda_events receipts
recorded in their params.xml (docs/pandda-campaign-design.md, section 10.2),
so existing campaigns are filled by ``manage.py backfill_campaign_events``,
not by a data migration here: the command is also the rebuild path after a
restore, and it reads files a migration has no business opening.
"""

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ccp4i2', '0026_campaign_uuids'),
    ]

    operations = [
        migrations.AddField(
            model_name='campaignsite',
            name='radius',
            field=models.FloatField(default=8.0),
        ),
        migrations.CreateModel(
            name='CampaignEvent',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('run_job_uuid', models.UUIDField(blank=True, null=True)),
                ('dtag', models.CharField(blank=True, default='', max_length=255)),
                ('event_idx', models.IntegerField()),
                ('site_idx', models.IntegerField(blank=True, null=True)),
                ('centroid_x', models.FloatField(blank=True, null=True)),
                ('centroid_y', models.FloatField(blank=True, null=True)),
                ('centroid_z', models.FloatField(blank=True, null=True)),
                ('hit_probability', models.FloatField(blank=True, null=True)),
                ('has_pose', models.BooleanField(default=False)),
                ('has_map', models.BooleanField(default=False)),
                ('cell_a', models.FloatField(blank=True, null=True)),
                ('cell_b', models.FloatField(blank=True, null=True)),
                ('cell_c', models.FloatField(blank=True, null=True)),
                ('cell_alpha', models.FloatField(blank=True, null=True)),
                ('cell_beta', models.FloatField(blank=True, null=True)),
                ('cell_gamma', models.FloatField(blank=True, null=True)),
                ('group', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='events', to='ccp4i2.projectgroup')),
                ('project', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='campaign_events', to='ccp4i2.project')),
                ('receipt', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='campaign_events', to='ccp4i2.job')),
            ],
            options={
                'constraints': [models.UniqueConstraint(fields=('receipt', 'site_idx', 'event_idx'), name='campaign_event_per_receipt')],
            },
        ),
    ]
