# Deploying API documentation using MKDocs + mkdocstrings

## Installation
```commandline
pip install mkdocs mkdocs-material mkdocstrings[python] mkdocstrings pyyaml
```

## Deploying
Make sure the package is installed locally

```pip install -e src```

Generate API level documentation
```commandline
# generate docs
python docs/tools/automatic_api_generation.py --package PAMparametrizer
# (optional) also update mkdocs nav:
python docs/tools/automatic_api_generation.py --package PAMparametrizer --update-mkdocs
```

## Adaptating the defaults
Check `docs/tools/automatic_api_generation.py`