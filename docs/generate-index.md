{% set command = 'generateIndex' %}
# {{ command }}

::: moPepGen.cli.generate_index
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