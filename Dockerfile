# OncoUnify: loaders + reference web interface on Apache (Debian 12).
# Demo:        docker compose up --build        (http://localhost:8080/)
# Institution: see docs/INSTALL.md (Option A)
FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        apache2 libcgi-pm-perl libdbi-perl libdbd-sqlite3-perl \
        python3 python3-pandas python3-openpyxl sqlite3 \
 && rm -rf /var/lib/apt/lists/* \
 && a2enmod cgid

WORKDIR /opt/oncounify
COPY . /opt/oncounify

RUN install -d /var/www/oncounify/cgi-bin /var/www/oncounify/panel /data /reports /etc/oncounify \
 && install -m 0755 /opt/oncounify/*.cgi /var/www/oncounify/cgi-bin/ \
 && install -m 0644 /opt/oncounify/search.html /var/www/oncounify/panel/ \
 && install -m 0644 /opt/oncounify/docker/apache-oncounify.conf /etc/apache2/sites-available/000-default.conf \
 && chmod 0755 /opt/oncounify/docker/*.sh

ENV ONCOUNIFY_DB=/data/panels.db \
    ONCOUNIFY_REPORTS=/reports \
    ONCOUNIFY_DEMO=0 \
    ONCOUNIFY_LOAD=0

EXPOSE 80
VOLUME ["/data"]
ENTRYPOINT ["/opt/oncounify/docker/entrypoint.sh"]
CMD ["apache2ctl", "-D", "FOREGROUND"]
