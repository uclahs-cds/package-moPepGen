{% set command = 'parseFusionCatcher' %}
# {{ command }}

::: moPepGen.cli.parse_fusion_catcher
	handler: python
    selection:
      members: false
    rendering:
      show_root_heading: false
      show_source: false

{% include 'partials/_caution_on_reference_version.md' %}

```
{{ parse_vep_help(command) }}
```