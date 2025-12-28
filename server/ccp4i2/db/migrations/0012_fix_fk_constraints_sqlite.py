# Generated manually to fix FK constraints in SQLite
# SQLite doesn't support ALTER TABLE to change FK constraints, so we need to recreate tables
# This migration is a no-op on PostgreSQL and other databases

from django.db import connection, migrations


def noop(apps, schema_editor):
    """No operation - used for non-SQLite databases."""
    pass


def recreate_jobfloatvalue(apps, schema_editor):
    """Recreate JobFloatValue table with correct FK (SQLite only)."""
    if connection.vendor != 'sqlite':
        return  # Only run on SQLite

    with connection.cursor() as cursor:
        cursor.executescript("""
            CREATE TABLE "ccp4i2_jobfloatvalue_new" (
                "id" integer NOT NULL PRIMARY KEY AUTOINCREMENT,
                "value" real NOT NULL,
                "job_id" bigint NOT NULL REFERENCES "ccp4i2_job" ("id") DEFERRABLE INITIALLY DEFERRED,
                "key_id" varchar(50) NOT NULL REFERENCES "ccp4i2_jobvaluekey" ("name") DEFERRABLE INITIALLY DEFERRED
            );
            INSERT INTO "ccp4i2_jobfloatvalue_new" SELECT * FROM "ccp4i2_jobfloatvalue";
            DROP TABLE "ccp4i2_jobfloatvalue";
            ALTER TABLE "ccp4i2_jobfloatvalue_new" RENAME TO "ccp4i2_jobfloatvalue";
            CREATE UNIQUE INDEX "ccp4i2_jobfloatvalue_job_id_key_id" ON "ccp4i2_jobfloatvalue" ("job_id", "key_id");
            CREATE INDEX "ccp4i2_jobfloatvalue_job_id" ON "ccp4i2_jobfloatvalue" ("job_id");
            CREATE INDEX "ccp4i2_jobfloatvalue_key_id" ON "ccp4i2_jobfloatvalue" ("key_id");
        """)


def recreate_jobcharvalue(apps, schema_editor):
    """Recreate JobCharValue table with correct FK (SQLite only)."""
    if connection.vendor != 'sqlite':
        return  # Only run on SQLite

    with connection.cursor() as cursor:
        cursor.executescript("""
            CREATE TABLE "ccp4i2_jobcharvalue_new" (
                "id" integer NOT NULL PRIMARY KEY AUTOINCREMENT,
                "value" varchar(255) NOT NULL,
                "job_id" bigint NOT NULL REFERENCES "ccp4i2_job" ("id") DEFERRABLE INITIALLY DEFERRED,
                "key_id" varchar(50) NOT NULL REFERENCES "ccp4i2_jobvaluekey" ("name") DEFERRABLE INITIALLY DEFERRED
            );
            INSERT INTO "ccp4i2_jobcharvalue_new" SELECT * FROM "ccp4i2_jobcharvalue";
            DROP TABLE "ccp4i2_jobcharvalue";
            ALTER TABLE "ccp4i2_jobcharvalue_new" RENAME TO "ccp4i2_jobcharvalue";
            CREATE UNIQUE INDEX "ccp4i2_jobcharvalue_job_id_key_id" ON "ccp4i2_jobcharvalue" ("job_id", "key_id");
            CREATE INDEX "ccp4i2_jobcharvalue_job_id" ON "ccp4i2_jobcharvalue" ("job_id");
            CREATE INDEX "ccp4i2_jobcharvalue_key_id" ON "ccp4i2_jobcharvalue" ("key_id");
        """)


class Migration(migrations.Migration):

    dependencies = [
        ('ccp4i2', '0011_fix_jobcharvalue_key_fk'),
    ]

    operations = [
        migrations.RunPython(recreate_jobfloatvalue, noop),
        migrations.RunPython(recreate_jobcharvalue, noop),
    ]
