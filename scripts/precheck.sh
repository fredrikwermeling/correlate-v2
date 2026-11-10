#!/bin/sh
# Checks run before every push (also installed as .git/hooks/pre-push):
# syntax, the statistics against scipy references, and a lint pass that
# stops on real errors (undefined names, duplicate keys and the like).
# The linter is optional: npm install eslint@8.57.0 --prefix scripts/.lint
cd "$(dirname "$0")/.." || exit 1
node --check app.js || exit 1
node --check tour.js || exit 1
node scripts/stat_check.js > /dev/null || { node scripts/stat_check.js; exit 1; }
ESLINT=scripts/.lint/node_modules/eslint/bin/eslint.js
if [ -f "$ESLINT" ]; then
    node "$ESLINT" --no-eslintrc -c scripts/eslint_check.json --quiet app.js tour.js || exit 1
else
    echo "precheck: linter not installed, lint skipped"
fi
echo "precheck: ok"
