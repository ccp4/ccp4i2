# Generated manually to fix FK constraints in SQLite
# SQLite doesn't support ALTER TABLE to change FK constraints, so we need to recreate tables

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('ccp4x', '0011_fix_jobcharvalue_key_fk'),
    ]

    operations = [
        # Recreate JobFloatValue with correct FK reference
        migrations.RunSQL(
            sql="""
                -- Create new table with correct FK
                CREATE TABLE "ccp4x_jobfloatvalue_new" (
                    "id" integer NOT NULL PRIMARY KEY AUTOINCREMENT,
                    "value" real NOT NULL,
                    "job_id" bigint NOT NULL REFERENCES "ccp4x_job" ("id") DEFERRABLE INITIALLY DEFERRED,
                    "key_id" varchar(50) NOT NULL REFERENCES "ccp4x_jobvaluekey" ("name") DEFERRABLE INITIALLY DEFERRED
                );
                -- Copy data from old table
                INSERT INTO "ccp4x_jobfloatvalue_new" SELECT * FROM "ccp4x_jobfloatvalue";
                -- Drop old table
                DROP TABLE "ccp4x_jobfloatvalue";
                -- Rename new table
                ALTER TABLE "ccp4x_jobfloatvalue_new" RENAME TO "ccp4x_jobfloatvalue";
                -- Recreate indexes
                CREATE UNIQUE INDEX "ccp4x_jobfloatvalue_job_id_key_id" ON "ccp4x_jobfloatvalue" ("job_id", "key_id");
                CREATE INDEX "ccp4x_jobfloatvalue_job_id" ON "ccp4x_jobfloatvalue" ("job_id");
                CREATE INDEX "ccp4x_jobfloatvalue_key_id" ON "ccp4x_jobfloatvalue" ("key_id");
            """,
            reverse_sql="""
                -- Reverse: Create table with old FK
                CREATE TABLE "ccp4x_jobfloatvalue_new" (
                    "id" integer NOT NULL PRIMARY KEY AUTOINCREMENT,
                    "value" real NOT NULL,
                    "job_id" bigint NOT NULL REFERENCES "ccp4x_job" ("id") DEFERRABLE INITIALLY DEFERRED,
                    "key_id" varchar(50) NOT NULL REFERENCES "ccp4x_jobvaluekeys" ("name") DEFERRABLE INITIALLY DEFERRED
                );
                INSERT INTO "ccp4x_jobfloatvalue_new" SELECT * FROM "ccp4x_jobfloatvalue";
                DROP TABLE "ccp4x_jobfloatvalue";
                ALTER TABLE "ccp4x_jobfloatvalue_new" RENAME TO "ccp4x_jobfloatvalue";
                CREATE UNIQUE INDEX "ccp4x_jobfloatvalue_job_id_key_id" ON "ccp4x_jobfloatvalue" ("job_id", "key_id");
                CREATE INDEX "ccp4x_jobfloatvalue_job_id" ON "ccp4x_jobfloatvalue" ("job_id");
                CREATE INDEX "ccp4x_jobfloatvalue_key_id" ON "ccp4x_jobfloatvalue" ("key_id");
            """,
        ),
        # Recreate JobCharValue with correct FK reference
        migrations.RunSQL(
            sql="""
                -- Create new table with correct FK
                CREATE TABLE "ccp4x_jobcharvalue_new" (
                    "id" integer NOT NULL PRIMARY KEY AUTOINCREMENT,
                    "value" varchar(255) NOT NULL,
                    "job_id" bigint NOT NULL REFERENCES "ccp4x_job" ("id") DEFERRABLE INITIALLY DEFERRED,
                    "key_id" varchar(50) NOT NULL REFERENCES "ccp4x_jobvaluekey" ("name") DEFERRABLE INITIALLY DEFERRED
                );
                -- Copy data from old table
                INSERT INTO "ccp4x_jobcharvalue_new" SELECT * FROM "ccp4x_jobcharvalue";
                -- Drop old table
                DROP TABLE "ccp4x_jobcharvalue";
                -- Rename new table
                ALTER TABLE "ccp4x_jobcharvalue_new" RENAME TO "ccp4x_jobcharvalue";
                -- Recreate indexes
                CREATE UNIQUE INDEX "ccp4x_jobcharvalue_job_id_key_id" ON "ccp4x_jobcharvalue" ("job_id", "key_id");
                CREATE INDEX "ccp4x_jobcharvalue_job_id" ON "ccp4x_jobcharvalue" ("job_id");
                CREATE INDEX "ccp4x_jobcharvalue_key_id" ON "ccp4x_jobcharvalue" ("key_id");
            """,
            reverse_sql="""
                -- Reverse: Create table with old FK
                CREATE TABLE "ccp4x_jobcharvalue_new" (
                    "id" integer NOT NULL PRIMARY KEY AUTOINCREMENT,
                    "value" varchar(255) NOT NULL,
                    "job_id" bigint NOT NULL REFERENCES "ccp4x_job" ("id") DEFERRABLE INITIALLY DEFERRED,
                    "key_id" varchar(50) NOT NULL REFERENCES "ccp4x_jobvaluekeys" ("name") DEFERRABLE INITIALLY DEFERRED
                );
                INSERT INTO "ccp4x_jobcharvalue_new" SELECT * FROM "ccp4x_jobcharvalue";
                DROP TABLE "ccp4x_jobcharvalue";
                ALTER TABLE "ccp4x_jobcharvalue_new" RENAME TO "ccp4x_jobcharvalue";
                CREATE UNIQUE INDEX "ccp4x_jobcharvalue_job_id_key_id" ON "ccp4x_jobcharvalue" ("job_id", "key_id");
                CREATE INDEX "ccp4x_jobcharvalue_job_id" ON "ccp4x_jobcharvalue" ("job_id");
                CREATE INDEX "ccp4x_jobcharvalue_key_id" ON "ccp4x_jobcharvalue" ("key_id");
            """,
        ),
    ]
