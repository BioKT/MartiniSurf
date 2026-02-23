# MartiniSurf Docs (Sphinx + Read the Docs)

This folder is prepared for Read the Docs deployment.

## Local preview

1. Install docs dependencies:

```bash
pip install -r docs/requirements.txt
```

2. Build HTML docs:

```bash
make -C docs html
```

3. Open locally:

- `docs/_build/html/index.html`

## Read the Docs setup

Configuration is already in:

- `.readthedocs.yaml`
- `docs/conf.py`
- `docs/requirements.txt`

When your GitHub repository becomes public:

1. Go to https://readthedocs.org/
2. Import your project from GitHub.
3. Select your repository and branch.
4. Build docs.

RTD will automatically use `.readthedocs.yaml` and publish the site.

## Notes

- Logo used in docs: `logo.png` from repository root.
- Main docs entrypoint: `docs/index.md`.
- Full user guide page: `docs/USER_GUIDE.md`.
