#!/usr/bin/perl
use strict;
use warnings;
use utf8;

binmode(STDOUT, ':encoding(UTF-8)');

print "Content-Type: text/html; charset=UTF-8\r\n";
print "Cache-Control: no-store, no-cache, must-revalidate, max-age=0\r\n";
print "Pragma: no-cache\r\n\r\n";

print <<'HTML';
<!DOCTYPE html>
<html lang="ja">
<head>
  <meta charset="UTF-8">
  <title>Logout</title>
  <style>
    body {
      margin: 0;
      padding: 2rem;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f3f4f6;
      color: #111827;
    }
    .card {
      max-width: 720px;
      margin: 0 auto;
      background: #fff;
      border-radius: 16px;
      padding: 1.6rem 1.8rem;
      box-shadow: 0 12px 30px rgba(15, 23, 42, 0.12);
    }
    h1 { margin-top: 0; }
    .note { color: #6b7280; line-height: 1.7; }
    .btn {
      display: inline-block;
      margin-top: 1rem;
      padding: 0.55rem 1rem;
      border-radius: 999px;
      text-decoration: none;
      background: linear-gradient(135deg, #0ea5e9, #2563eb);
      color: #fff;
      font-weight: 600;
    }
  </style>
</head>
<body>
  <div class="card">
    <h1>Logout</h1>
    <p class="note">
      With HTTP Basic authentication, browsers cache credentials differently.<br>
      This page attempts to invalidate the cached credentials, but in some browsers the credentials may persist until the tab or the browser is closed.
    </p>
    <p class="note" id="status">Signing out&hellip;</p>
    <a class="btn" href="/panel/search.html">Return to search</a>
  </div>

  <script>
    (async function() {
      const status = document.getElementById('status');
      try {
        if (document.execCommand) {
          try { document.execCommand('ClearAuthenticationCache'); } catch (e) {}
        }
        try {
          await fetch('/panel/search.html', {
            method: 'GET',
            cache: 'no-store',
            headers: {
              'Authorization': 'Basic ' + btoa('logout:logout')
            }
          });
        } catch (e) {}
        status.textContent = 'Sign-out attempted. If the authentication dialog appears on the next request, you have been logged out.';
      } catch (e) {
        status.textContent = 'Sign-out attempted. If your browser still remembers the credentials, close this tab or the browser.';
      }
    })();
  </script>
</body>
</html>
HTML
