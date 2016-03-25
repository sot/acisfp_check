===========================
FP_TEMP temperatures check
===========================
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
Date start            {{proc.datestart}}
Date stop             {{proc.datestop}}
FP_TEMP status        {%if red_viols.fptemp or fp_sens_viols.fptemp%}:red:`NOT OK`{% else %}OK{% endif%} 
{% if opt.loaddir %}
Load directory        {{opt.loaddir}}
{% endif %}
Run time              {{proc.run_time}} by {{proc.run_user}}
Run log               `<run.dat>`_
Temperatures          `<temperatures.dat>`_
States                `<states.dat>`_
====================  =============================================

{% if red_viols.fptemp  %}
FP_TEMP -114 deg C Violations
------------------------------
=====================  =====================  ==================  ==================
Date start             Date stop              Max temperature     Obsids
=====================  =====================  ==================  ==================
{% for viol in red_viols.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{viol.maxtemp|floatformat:2}}            {{viol.obsid}}
{% endfor %}
=====================  =====================  ==================  ==================
{% else %}
No -114 deg C FP_TEMP Violations
{% endif %}


{% if fp_sens_viols.fptemp  %}
FP_TEMP -118.7 deg C Violations for FP TEMP Sensitive Observations
-------------------------------------------------------------------
=====================  =====================  ==================  ==================
Date start             Date stop              Max temperature     OBSID
=====================  =====================  ==================  ==================
{% for viol in fp_sens_viols.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{viol.maxtemp|floatformat:2}}             {{viol.obsid}}
{% endfor %}
=====================  =====================  ==================  ==================
{% else %}
No Focal Plane Sensitive Observation -118.7 deg C FP_TEMP Violations
{% endif %}



{% if cti_viols.fptemp  %}
FP_TEMP -118.7 deg C Violations for Perigee Passages
-------------------------------------------------------------------
=====================  =====================  ==================
Date start             Date stop              Max temperature
=====================  =====================  ==================
{% for viol in cti_viols.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{viol.maxtemp|floatformat:2}}
{% endfor %}
=====================  =====================  ==================
{% else %}
No ECS Observation -118.7 deg C FP_TEMP Violations
{% endif %}



.. image:: {{plots.fptemp.filename}}
.. image:: {{plots.pow_sim.filename}}

=========================
FP_TEMP Model Validation
=========================

MSID quantiles
---------------

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
   {{viol.msid}},{{viol.quant}},{{viol.value}},{{viol.limit|floatformat:2}}
{% endfor%}

{% else %}
No Validation Violations
{% endif %}


{% for plot in plots_validation %}
{{ plot.msid }}
-----------------------


Red = telemetry, blue = model

.. image:: {{plot.lines}}
.. image:: {{plot.histlog}}
.. image:: {{plot.histlin}}

{% endfor %}

ADDITIONAL PLOTS
-----------------------

Additional plots of FPTEMP vs TIME for different temerature ranges

.. image:: fptempM120toM119.png
.. image:: fptempM120toM90.png
