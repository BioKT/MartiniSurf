# MartiniSurf Protein Designer

Streamlit interface for the MartiniSurf protein workflow. It generates design/setup files only; molecular dynamics execution is intentionally left outside the app.

## Local Run

```bash
pip install -r requirements-streamlit.txt
pip install -e .
streamlit run streamlit_app.py
```

The app writes temporary runs to `streamlit_runs/`, shows the generated command, captures logs, and creates a downloadable `Simulation_Files.zip`.
When launched from the repository, the app automatically adds the active Python environment and local `.venv/bin` to `PATH`, so tools installed there, such as `martinize2`, are available to the MartiniSurf runner.

## Streamlit Cloud

Use `streamlit_app.py` as the entrypoint and `requirements-streamlit.txt` as the dependency file. External tools such as `martinize2` and GROMACS must be available in the deployed environment for the corresponding workflow options to run.

The molecular viewer is embedded in the results panel and uses 3Dmol.js in the browser to inspect generated `.pdb` and `.gro` files.
