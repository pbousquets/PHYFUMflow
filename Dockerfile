FROM pbousquets/phyfum_base

RUN pip install phyfum

ENTRYPOINT ["phyfum"]
