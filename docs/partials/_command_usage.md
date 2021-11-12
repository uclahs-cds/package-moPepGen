{% for action in actions %}
**{{ ', '.join(action.option_strings) }}**{% if action.metavar != none %} ***{{action.metavar|e}}***{% endif %}{% if action.type != none and action.type.__name__ != "bool" %} *{{action.type.__name__}}*{% endif %}<br/>
<p style="margin-left: 15px;">
{{ action.help }}
{% if action.default != none and action.type.__name__ != "bool" and action.default != "==SUPPRESS==" %} {{ action.type.__name__ }} <br/>Default: {{ action.default }} {% endif %}
{% if action.choices != none %} <br/>Choices: {{ action.choices }} {% endif %}
</p>
{% endfor %}