# BRC Convenience Webapps
This repository contains marimo-made interactive web apps to process BRC-specific data. The sections below detail how to add new
apps.


## Including data or assets

To include data or assets in your notebooks, add them to the `public/` directory.

For example, the `apps/charts.py` notebook loads an image asset from the `public/` directory.

```markdown
<img src="public/logo.png" width="200" />
```

And the `notebooks/penguins.py` notebook loads a CSV dataset from the `public/` directory.

```python
import polars as pl
df = pl.read_csv(mo.notebook_location() / "public" / "penguins.csv")
```

## 🎨 Templates

This repository includes several templates for the generated site:

1. `index.html.j2` (default): A template with styling and a footer
2. `bare.html.j2`: A minimal template with basic styling
3. `tailwind.html.j2`: A minimal and lean template using Tailwind CSS

To use a specific template, pass the `--template` parameter to the build script:

```bash
uv run .github/scripts/build.py --template templates/tailwind.html.j2
```

You can also create your own custom templates. See the [templates/README.md](templates/README.md) for more information.

## 🧪 Testing

To test the export process, run `.github/scripts/build.py` from the root directory.

```bash
uv run .github/scripts/build.py
```

This will export all notebooks in a folder called `_site/` in the root directory. Then to serve the site, run:

```bash
python -m http.server -d _site
```

This will serve the site at `http://localhost:8000`.
