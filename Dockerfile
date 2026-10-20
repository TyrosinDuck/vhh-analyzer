# VHH Analyzer — container image.
#
# Based on micromamba so that HMMER (the ANARCI backend) can be installed from
# bioconda; a pip-only image would fall back to the regex annotator.
# Build:  docker build -t vhh-analyzer .
# Run:    docker run -p 8501:8501 vhh-analyzer

FROM mambaorg/micromamba:1.5-jammy

USER root
RUN useradd -m -u 1000 appuser || true
USER $MAMBA_USER

# HMMER + ANARCI come from bioconda; everything else from pip
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 hmmer anarci && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt /tmp/requirements.txt
RUN micromamba run -n base pip install --no-cache-dir -r /tmp/requirements.txt

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER vhh_core.py app.py /app/

EXPOSE 8501
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s \
    CMD curl -fsS http://localhost:8501/_stcore/health || exit 1

ENV STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

ENTRYPOINT ["micromamba", "run", "-n", "base", "streamlit", "run", "app.py"]
