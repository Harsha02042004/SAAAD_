import sqlite3
import re

conn = sqlite3.connect("data/sialic.db")   # ⚠️ adjust path if needed
cursor = conn.cursor()

def slugify(name):
    return re.sub(r'[^a-z0-9]+', '-', name.lower()).strip('-')

cursor.execute("SELECT compound_id, entry_name FROM compounds")

for cid, name in cursor.fetchall():
    slug = slugify(name)
    cursor.execute(
        "UPDATE compounds SET slug=? WHERE compound_id=?",
        (slug, cid)
    )

conn.commit()
conn.close()

print("✅ Slugs generated successfully")