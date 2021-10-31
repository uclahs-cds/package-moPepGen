{% set command = 'splitDatabase' %}
# {{ command }}

::: moPepGen.cli.call_variant_peptide
	handler: python
    selection:
      members: false
    rendering:
      show_root_heading: false
      show_source: false

```
{{ parse_vep_help(command) }}
```