{% set command = 'callVariant' %}
# {{ command }}

::: moPepGen.cli.call_variant_peptide
	handler: python
    selection:
      members: false
    rendering:
      show_root_heading: false
      show_source: false

## Arguments

{% with actions=get_arg_data(command) %}
{% include 'partials/_command_usage.md' %}
{% endwith %}