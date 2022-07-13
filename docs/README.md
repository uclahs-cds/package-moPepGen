## MKDocs

Documentations are written in markdown and converted to html by [MKDocs](https://www.mkdocs.org/).

## Adding Pages

To add a new documentation page, first create a '.md' file under this directory. Next go to the 'mkdocs.yml' file under the root directory of the repo, and add a key-value pair under `nav`, and after rendering, a link to the new added page will show up in the navbar.

## Render locally

The online documentation is rendered automatically by github actions, it is still useful to see the changes in real time. This can be done by running `mkdocs` locally. moPepGen, mkdocs and several dependencies need to be installed first.

```bash
conda env create --name moPepGen python=3.8.11
conda activate moPepGen
pip install mkdocs==1.2.3 jinja2==3.0.0 mkdocstrings==0.16.2 mkdocs-macros-plugin==0.6.0
pip install .
```

Use the command below to start a development server and open http://127.0.0.1:8000/ to see the change in real time.

```bash
mkdocs serve
```
