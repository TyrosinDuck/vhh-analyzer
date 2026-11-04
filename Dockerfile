# VHH Analyzer, container image.
#
# Based on micromamba so that HMMER (the ANARCI backend) can be installed from
# bioconda; a pip-only image would fall back to the regex annotator.
# Build:  docker build -t vhh-analyzer .
# Run:    docker run -p 8501:8501 vhh-analyzer

FROM mambaorg/micromamba:1.5-jammy

USER root
RUN useradd -m -u 1000 appuser || true
USER $MAMBA_USER

# Python and HMMER come from conda-forge and bioconda; everything else,
# ANARCI included, from pip. The Python version is the one requirements.txt
# was compiled for (see its header): the lockfile pins packages that do not
# install on older interpreters, so the two have to move together, and
# tests/test_metadata.py fails if they do not. ANARCI used to come from
# bioconda as well, which put a second, unpinned copy beside the one the lock
# names.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.14 hmmer && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt /tmp/requirements.txt
RUN micromamba run -n base pip install --no-cache-dir -r /tmp/requirements.txt

WORKDIR /app
# The whole project, filtered by .dockerignore, rather than a list of files.
# The list this replaces named two files while app.py imports a dozen modules
# and reads sections/, data/, MANUAL.md and .streamlit/; nobody updated it as
# modules were added, and the image failed on its first import. What stays
# out is decided in .dockerignore, which only has to change when something
# should be excluded, not every time something is added.
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app/

EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s \
    CMD curl -fsS http://localhost:8501/_stcore/health || exit 1

ENV STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

ENTRYPOINT ["micromamba", "run", "-n", "base", "streamlit", "run", "app.py"]
