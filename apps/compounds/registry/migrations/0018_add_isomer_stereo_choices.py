"""
Add isomer_1 through isomer_20 to stereo_comment choices.

For compounds with multiple stereo centers where the exact configuration
is not known at registration time.
"""

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('registry', '0017_populate_inchi'),
    ]

    operations = [
        migrations.AlterField(
            model_name='compound',
            name='stereo_comment',
            field=models.CharField(
                choices=[
                    ('unset', 'Unset'),
                    ('achiral', 'Achiral'),
                    ('racemic', 'Racemic mixture'),
                    ('single_unknown', 'Single enantiomer, configuration unknown'),
                    ('r_enantiomer', 'R enantiomer'),
                    ('s_enantiomer', 'S enantiomer'),
                    ('non_racemic_mixture', 'Non-racemic stereoisomer mixture'),
                    ('four_diastereomers', 'Mixture of 4 diastereoisomers'),
                    ('two_diastereomers', 'Mixture of 2 diastereoisomers'),
                    ('single_diastereomer_unknown', 'Single diastereoisomer, configuration unknown'),
                    ('rr_diastereomer', 'RR diastereoisomer'),
                    ('rs_diastereomer', 'RS diastereoisomer'),
                    ('sr_diastereomer', 'SR diastereoisomer'),
                    ('ss_diastereomer', 'SS diastereoisomer'),
                    ('epimer_mixture', 'Mixture of epimers'),
                    ('ez_mixture', 'Mixture of E and Z isomers'),
                    ('e_isomer', 'E isomer'),
                    ('z_isomer', 'Z isomer'),
                    ('isomer_1', 'Isomer 1'),
                    ('isomer_2', 'Isomer 2'),
                    ('isomer_3', 'Isomer 3'),
                    ('isomer_4', 'Isomer 4'),
                    ('isomer_5', 'Isomer 5'),
                    ('isomer_6', 'Isomer 6'),
                    ('isomer_7', 'Isomer 7'),
                    ('isomer_8', 'Isomer 8'),
                    ('isomer_9', 'Isomer 9'),
                    ('isomer_10', 'Isomer 10'),
                    ('isomer_11', 'Isomer 11'),
                    ('isomer_12', 'Isomer 12'),
                    ('isomer_13', 'Isomer 13'),
                    ('isomer_14', 'Isomer 14'),
                    ('isomer_15', 'Isomer 15'),
                    ('isomer_16', 'Isomer 16'),
                    ('isomer_17', 'Isomer 17'),
                    ('isomer_18', 'Isomer 18'),
                    ('isomer_19', 'Isomer 19'),
                    ('isomer_20', 'Isomer 20'),
                ],
                default='unset',
                max_length=70,
            ),
        ),
    ]
