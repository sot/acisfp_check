==========================
FP_TEMP temperature check
==========================
.. role:: red


Summary
--------         
.. class:: borderless

====================  =============================================
Date start            2012:142:08:21:00.000
Date stop             2012:149:12:04:37.306
1DPAMZT status        OK (limit =  C)
Run time              Wed Jun  6 15:51:35 2012 by gregg
Run log               `<run.dat>`_
Temperatures          `<temperatures.dat>`_
States                `<states.dat>`_
====================  =============================================

No 1DPAMZT Violations

.. image:: 
.. image:: pow_sim.png

=======================
DPA Model Validation
=======================

MSID quantiles
---------------

Note: DPA quantiles are calculated using only points where 1DPAMZT > 20 degC.

.. csv-table:: 
   :header: "MSID", "1%", "5%", "16%", "50%", "84%", "95%", "99%"
   :widths: 15, 10, 10, 10, 10, 10, 10, 10

   PITCH,-1.205,0.089,0.184,0.256,0.317,0.418,1.315
   TSCPOS,-1,-1,-1,-1,0,0,0


No Validation Violations


PITCH
-----------------------

Note: DPA residual histograms include only points where 1DPAMZT > 20 degC.

Red = telemetry, blue = model

.. image:: pitch_valid.png
.. image:: pitch_valid_hist_log.png
.. image:: pitch_valid_hist_lin.png

TSCPOS
-----------------------

Note: DPA residual histograms include only points where 1DPAMZT > 20 degC.

Red = telemetry, blue = model

.. image:: tscpos_valid.png
.. image:: tscpos_valid_hist_log.png
.. image:: tscpos_valid_hist_lin.png


ADDITIONAL PLOTS
-----------------------

Note: DPA residual histograms include only points where 1DPAMZT > 20 degC.

Additional plots of FPTEMP vs TIME for different temerature ranges

.. image:: fptempM120toM119.png
.. image:: fptempM120toM90.png
