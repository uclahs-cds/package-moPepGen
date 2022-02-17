{% set command = 'decoyFasta' %}
# {{ command }}

::: moPepGen.cli.decoy_fasta
	handler: python
    selection:
      members: false
    rendering:
      show_root_heading: false
      show_source: false

## Usage

```
{{ get_arg_usage(command) }}
```

## Arguments

{% with actions=get_arg_data(command) %}
{% include 'partials/_command_usage.md' %}
{% endwith %}
