===========================
FP_TEMP temperatures check
===========================
.. role:: red

{% if proc.errors %}
Processing Errors
-----------------
.. class:: red
{% endif %}

{% if bsdir %}

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
Earth Solid Angles    `<earth_solid_angles.dat>`_
States                `<states.dat>`_
====================  =============================================

{% if viols.ACIS_I.fptemp %}
ACIS-I FP_TEMP -112 deg C Violations
------------------------------------
=====================  =====================  ==================  ==================
Date start             Date stop              Max temperature     Obsids
=====================  =====================  ==================  ==================
{% for viol in viols.ACIS_I.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{"%.2f"|format(viol.maxtemp)}}             {{viol.obsid}}
{% endfor %}
=====================  =====================  ==================  ==================
{% else %}
No ACIS-I -112 deg C FP_TEMP Violations
{% endif %}


{% if viols.ACIS_S.fptemp %}
ACIS-S FP_TEMP -111 deg C Violations
------------------------------------
=====================  =====================  ==================  ==================
Date start             Date stop              Max temperature     Obsids
=====================  =====================  ==================  ==================
{% for viol in viols.ACIS_S.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{"%.2f"|format(viol.maxtemp)}}             {{viol.obsid}}
{% endfor %}
=====================  =====================  ==================  ==================
{% else %}
No ACIS-S -111 deg C FP_TEMP Violations
{% endif %}

{% if viols.ecs.fptemp %}
Science Orbit ECS -119.5 deg C Violations
-----------------------------------------
=====================  =====================  ==================  ==================
Date start             Date stop              Max temperature     OBSID
=====================  =====================  ==================  ==================
{% for viol in viols.ecs.fptemp %}
{{viol.datestart}}  {{viol.datestop}}  {{"%.2f"|format(viol.maxtemp)}}             {{viol.obsid}}
{% endfor %}
=====================  =====================  ==================  ==================
{% else %}
No Science Orbit ECS -119.5 deg C FP_TEMP Violations
{% endif %}

.. image:: {{plots.acisfp_3.filename}}
.. image:: {{plots.pow_sim.filename}}
.. image:: {{plots.roll_taco.filename}}

{% endif %}

{% if not pred_only %}

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

{% if plot.msid == "ccd_count" %}

CCD/FEP Count
-------------

.. image:: {{plot.lines}}

{% elif plot.msid == "earthheat__fptemp" %}

Earth Solid Angle
-----------------

.. image:: {{plot.lines}}

{% else %}

{{ plot.msid }}
-----------------------


Red = telemetry, blue = model

.. image:: {{plot.lines}}

Data for FPTEMP residual plots limited between -120.0 and -112.0 deg. C
-----------------------------------------------------------------------

.. image:: {{plot.hist}}

{% endif %}

{% endfor %}

{% endif %}

{% if bsdir %}

ADDITIONAL PLOTS
-----------------------

Additional plots of FPTEMP vs TIME for different temerature ranges

.. image:: fptempM120toM119.png
.. image:: fptempM120toM90.png

{% endif %}
