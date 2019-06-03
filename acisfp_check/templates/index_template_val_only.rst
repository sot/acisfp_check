=======================
FP_Temp temperatures check
=======================
.. role:: red

{% if proc.errors %}
Processing Errors
-----------------
.. class:: red
{% endif %}

Summary
--------         
.. class:: borderless

====================  =============================================
{% if bsdir %}
Load directory        {{bsdir}}
{% endif %}
Run time              {{proc.run_time}} by {{proc.run_user}}
Run log               `<run.dat>`_
====================  =============================================

=======================
FP Model Validation
=======================

MSID quantiles
---------------

Note: DPA quantiles are calculated using only points where 1DPAMZT > 20 degC.

.. csv-table:: 
   :header: "MSID", "1%", "5%", "16%", "50%", "84%", "95%", "99%"
   :widths: 15, 10, 10, 10, 10, 10, 10, 10

{% for plot in plots_validation %}
{% if plot.quant01 %}
   {{plot.msid}},{{plot.quant01}},{{plot.quant05}},{{plot.quant16}},{{plot.quant50}},{{plot.quant84}},{{plot.quant95}},{{plot.quant99}}
{% endif %}
{% endfor%}


{% if valid_viols %}
Validation Violations
---------------------

.. csv-table:: 
   :header: "MSID", "Quantile", "Value", "Limit"
   :widths: 15, 10, 10, 10

{% for viol in valid_viols %}
   {{viol.msid}},{{viol.quant}},{{viol.value}},{{"%.2f"|format(viol.limit)}}
{% endfor%}

{% else %}
No Validation Violations
{% endif %}
   
{% for plot in plots_validation %}
{{ plot.msid }}
-----------------------


Red = telemetry, blue = model

.. image:: {{plot.lines}}

Data for FPTEMP residual plots limited between -120.0 and -112.0 deg. C
-----------------------------------------------------------------------

.. image:: {{plot.hist}}

{% endfor %}
