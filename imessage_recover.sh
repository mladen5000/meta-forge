#!/usr/bin/env zsh
set -euo pipefail

# --------- CONFIG: adjust your window here ---------
START="2025-07-14 06:26:16"
END="2025-07-14 06:30:16"
KEYWORDS=("password" "otp" "code")
# ---------------------------------------------------

MSG_DIR="$HOME/Library/Messages"
WORK="$HOME/Documents"
DB="$WORK/chat_copy.db"
WAL="$WORK/chat_copy.db-wal"
SHM="$WORK/chat_copy.db-shm"
CSV="$WORK/imessage_window_export.csv"
CARVER="$WORK/carve_imessage_texts.py"
CAND="$WORK/carved_candidates.txt"

echo "[*] Copying DB/WAL/SHM"
cp -v "$MSG_DIR/chat.db" "$DB"
cp -v "$MSG_DIR/chat.db-wal" "$WAL" || true
cp -v "$MSG_DIR/chat.db-shm" "$SHM" || true

echo "[*] Exporting joined window CSV: $CSV"
sqlite3 "$DB" <<SQL
.mode csv
.headers on
WITH window AS (
  SELECT datetime('$START') AS start_ts, datetime('$END') AS end_ts
)
SELECT
  datetime(m.date/1000000000 + strftime('%s','2001-01-01'), 'unixepoch', 'localtime') as msg_date_local,
  CASE m.is_from_me WHEN 1 THEN 'me' ELSE 'other' END as direction,
  COALESCE(h.id,'') as handle,
  COALESCE(c.display_name,'') as chat_display_name,
  COALESCE(c.guid,'') as chat_guid,
  COALESCE(m.text,'') as text,
  m.ROWID as message_rowid
FROM message m
LEFT JOIN handle h ON h.ROWID=m.handle_id
LEFT JOIN chat_message_join cmj ON cmj.message_id=m.ROWID
LEFT JOIN chat c ON c.ROWID=cmj.chat_id
WHERE datetime(m.date/1000000000 + strftime('%s','2001-01-01'), 'unixepoch', 'localtime')
      BETWEEN (SELECT start_ts FROM window) AND (SELECT end_ts FROM window)
ORDER BY msg_date_local, m.ROWID;
SQL
sqlite3 "$DB" ".once $CSV" ".read /dev/stdin" < /dev/null || true

echo "[*] WAL strings dump (if any)"
if [[ -f "$WAL" ]]; then
  LC_ALL=C strings -a "$WAL" > "$WORK/wal_strings.txt"
  echo "  -> $WORK/wal_strings.txt"
fi

if [[ ! -f "$CARVER" ]]; then
  cat > "$CARVER" <<'PY'
# (paste the exact Python carver from section 2B here)
PY
  echo "[*] WROTE $CARVER — paste the carver code into it and rerun if empty"
fi

echo "[*] Running carver"
python3 "$CARVER" --db "$DB" --wal "$WAL" --max 500 --minlen 8 --out "$CAND" \
  --keywords "${KEYWORDS[@]}" || python3 "$CARVER" --db "$DB" --max 500 --minlen 8 --out "$CAND"

echo
echo "Done."
echo "  CSV:     $CSV"
echo "  CARVED:  $CAND"
[[ -f "$WORK/wal_strings.txt" ]] && echo "  WAL TXT: $WORK/wal_strings.txt"

