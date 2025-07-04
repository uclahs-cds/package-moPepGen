{% set command = 'callNovelORF' %}
# {{ command }}

::: moPepGen.cli.call_novel_orf
	handler: python
    selection:
      members: false
    rendering:
      show_root_heading: false
      show_source: false

{% include 'partials/_caution_on_reference_version.md' %}

{% include 'partials/_args_reference.md' %}

{% include 'partials/_args_codon_table.md' %}

## Usage

```
{{ get_arg_usage(command) }}
```

## Arguments

{% with actions=get_arg_data(command) %}
{% include 'partials/_command_usage.md' %}
{% endwith %}
