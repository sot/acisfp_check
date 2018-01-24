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
FP_TEMP status        {%if viols.ACIS_I.fptemp or viols.ACIS_S.fptemp %}:red:`NOT OK`{% else %}OK{% endif%} 
{% if bsdir %}
Load directory        {{bsdir}}
{% endif %}
Run time              {{proc.run_time}} by {{proc.run_user}}
Run log               `<run.dat>`_
Temperatures          `<temperatures.dat>`_
States                `<states.dat>`_
====================  =============================================

{% if ACIS_I_viols.fptemp %}
ACIS-I FP_TEMP -114 deg C Violations
------------------------------------
=====================  =====================  ==================  ==================
Date start             Date stop              Max temperature     Obsids
=====================  =====================  ==================  ==================
{% for viol in viols.ACIS_I.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{"%.2f"|format(viol.maxtemp)}}            {{viol.obsid}}
{% endfor %}
=====================  =====================  ==================  ==================
{% else %}
No ACIS-I -114 deg C FP_TEMP Violations
{% endif %}


{% if ACIS_S_viols.fptemp %}
ACIS-S FP_TEMP -112 deg C Violations
------------------------------------
=====================  =====================  ==================  ==================
Date start             Date stop              Max temperature     Obsids
=====================  =====================  ==================  ==================
{% for viol in viols.ACIS_S.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{"%.2f"|format(viol.maxtemp)}}            {{viol.obsid}}
{% endfor %}
=====================  =====================  ==================  ==================
{% else %}
No ACIS-S -112 deg C FP_TEMP Violations
{% endif %}


{% if viols.fp_sens.fptemp %}
FP TEMP Sensitive, -118.7 deg. C Preference Not Met:
-------------------------------------------------------------------
=====================  =====================  ==================  ==================
Date start             Date stop              Max temperature     OBSID
=====================  =====================  ==================  ==================
{% for viol in viols.fp_sens.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{"%.2f"|format(viol.maxtemp)}}             {{viol.obsid}}
{% endfor %}
=====================  =====================  ==================  ==================
{% else %}
No Focal Plane Sensitive Observation -118.7 deg C FP_TEMP Preferences Unmet
{% endif %}



{% if viols.cti.fptemp %}
FP_TEMP -118.7 deg C Violations for Perigee Passages
-------------------------------------------------------------------
=====================  =====================  ==================
Date start             Date stop              Max temperature
=====================  =====================  ==================
{% for viol in viols.cti.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{"%.2f"|format(viol.maxtemp)}}
{% endfor %}
=====================  =====================  ==================
{% else %}
No ECS Observation -118.7 deg C FP_TEMP Violations
{% endif %}



.. image:: {{plots.acisfp_3.filename}}
.. image:: {{plots.pow_sim.filename}}
.. image:: {{plots.roll.filename}}

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
{% endfor %}

{% if valid_viols %}
Validation Violations
---------------------

.. csv-table:: 
   :header: "MSID", "Quantile", "Value", "Limit"
   :widths: 15, 10, 10, 10


{% for viol in valid_viols %}
   {{viol.msid}},{{viol.quant}},{{viol.value}},{{"%.2f"|format(viol.limit)}}
{% endfor %}

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

.. image:: {{plot.histlog}}
.. image:: {{plot.histlin}}

{% endfor %}

ADDITIONAL PLOTS
-----------------------

Additional plots of FPTEMP vs TIME for different temerature ranges

.. image:: fptempM120toM119.png
.. image:: fptempM120toM90.png
