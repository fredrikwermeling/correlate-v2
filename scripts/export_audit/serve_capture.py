#!/usr/bin/env python3
"""Serve the V2 repo and accept POST /save?name= chunks (export-audit capture).

Automated Chrome cannot complete downloads, so the audit shims in the page
POST exported blobs here instead; each POST APPENDS to out/<name>, so use one
unique name per capture.
"""
import http.server, os, urllib.parse

ROOT = "/Users/fredrikwermeling/Documents/correlate_v2"
OUT = os.path.join(ROOT, "scripts", "export_audit", "out")

class H(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=ROOT, **k)
    def do_POST(self):
        q = urllib.parse.urlparse(self.path)
        if q.path != "/save":
            self.send_error(404); return
        name = urllib.parse.parse_qs(q.query).get("name", ["capture"])[0]
        name = os.path.basename(name)
        os.makedirs(OUT, exist_ok=True)
        n = int(self.headers.get("Content-Length", 0))
        data = self.rfile.read(n)
        with open(os.path.join(OUT, name), "ab") as f:
            f.write(data)
        self.send_response(200); self.end_headers(); self.wfile.write(b"ok")
    def log_message(self, *a): pass

http.server.ThreadingHTTPServer(("127.0.0.1", 8642), H).serve_forever()
