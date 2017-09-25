===============================================================
**How to calculate how light scatters from a Source to a Site**
===============================================================

**Summary and References**
------------------------------------
Our model of light propagation is based on three papers.

(1) Falchi, F., P. Cinzano, D. Duriscoe, C.C.M. Kyba, C.D. Elvidge, K. Baugh, B.A. Portnov, N.A. Rybnikova and R. Furgoni, 2016. The new world atlas of artificial night sky brightness. Sci. Adv. 2.
(2) Cinzano, P., F. Falchi, C.D. Elvidge and  K.E. Baugh, 2000. The artificial night sky brightness mapped from DMSP satellite Operational Linescan System measurements. Mon. Not. R. Astron. Soc. 318.
(3) Garstang, R.H., 1989. Night-sky brightness at observatories and sites. Pub. Astron. Soc. Pac. 101.

These papers are the foundation for the part of our code that calculates how light scatters from a source (C) to a site (O for observer). In our case, values of light emission from the source are sensed by VIIRS and recorded as pixel values, and each pixel area over our study area is a site. This explanation will cover the inner workings of the function ``fsum_single()``, which takes as arguments parameters that govern how light scatters from a C to O, and outputs a ratio that determines the total amount of light that arrives at the line of sight of viewer standing at O.

The following adapted diagram from (2) shows the relavant geometry parameters that govern how light scatters from a C to a O:

.. image:: ../_static/cinzano_diagram.png

``fsum_single()`` **Code Organization**
---------------------------------------

Each equation in ``fsum_single()`` that is taken from any of the papers that detail our model for light pollution are referenced as follows:

.. code-block:: python
	
    # s, Distance from source to scattering CQ, REF 2, Appendix A (A1), p. 656
    # equation is wrong in Ref 2 (Cinzano). Changed Chi to theta according to Ref 3, p. 308, Equation 7 (Garstang)
    s_CQ = sqrt((u_OQ - l_OC)**2 + 4*u_OQ*l_OC*sin(theta/2)**2) # km

    # h, Height of scattering above Earth reference surface Q, REF 2, Appendix A (A1), p. 656
    h_Q = R_T*(sqrt(1 + (u_OQ**2 + 2*u_OQ*R_T*cos(zen))/R_T**2) - 1) # km

    # phi, Elevation angle of emission over line of sight from C to O (QCO), REF 2, Appendix A (A1), p. 656
    phi = arccos((l_OC**2 + s_CQ**2 - u_OQ**2)/(2*l_OC*s_CQ)) # radians
    phi_deg = phi*180/pi # degrees

Above each equation or constant, there is a comment with a short explanation of the variable, a reference, and page number. Starting out with the code, there are a lot of equations to understand and the task of digesting them all can seem daunting! never fear. You do not need to meticulously know the reason behind every equation unless you are interested in the physics/geometry of light scattering. It's just important to know where these equations come from and how to examine their meaning should you need to do so for debugging purposes. We will go over some of the most relevant parts of ``fsum_single()`` that you may need to change.


