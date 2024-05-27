package IBPSA "Library with models for building energy and control systems"
  extends Modelica.Icons.Package;

  package Fluid "Package with models for fluid flow systems"
    extends Modelica.Icons.Package;

    package Sensors "Package with sensor models"
      extends Modelica.Icons.SensorsPackage;

      model RelativeHumidity "Ideal one port relative humidity sensor"
        extends IBPSA.Fluid.Sensors.BaseClasses.PartialAbsoluteSensor;
        extends Modelica.Icons.RoundSensor;

        Modelica.Blocks.Interfaces.RealOutput phi(final unit="1", min=0)
          "Relative humidity in port medium"
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));

    protected
        Modelica.Units.SI.Temperature T "Temperature of the medium";
        Medium.MassFraction Xi[Medium.nXi](
          quantity=Medium.substanceNames[1:Medium.nXi]) "Mass fraction of the medium";
      equation
        Xi = inStream(port.Xi_outflow);
        T=Medium.temperature_phX(
            p=port.p,
            h=inStream(port.h_outflow),
            X=Xi);

        phi = IBPSA.Utilities.Psychrometrics.Functions.phi_pTX(
          p=port.p,
          T=T,
          X_w=Xi[1]);

      annotation (defaultComponentName="senRelHum",
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
              graphics={
              Line(points={{0,-70},{0,-100}}, color={0,0,127}),
              Text(
                extent={{-150,80},{150,120}},
                textString="%name",
                textColor={0,0,255}),
              Text(
                extent={{160,-30},{60,-60}},
                textColor={0,0,0},
                textString="phi"),
              Line(points={{70,0},{100,0}}, color={0,0,127}),
              Text(
                extent={{180,90},{60,40}},
                textColor={0,0,0},
                textString=DynamicSelect("", String(phi, leftJustified=false, significantDigits=2)))}),
        Documentation(info="<html>
<p>
This model outputs the relative humidity of the fluid connected to its port.
The sensor is ideal, i.e. it does not influence the fluid.
</p>
<p>
Note that this sensor can only be used with media that contain the variable <code>phi</code>,
which is typically the case for moist air models.
</p>
<p>
To measure relative humidity in a duct or pipe, use
<a href=\"modelica://IBPSA.Fluid.Sensors.RelativeHumidityTwoPort\">IBPSA.Fluid.Sensors.RelativeHumidityTwoPort</a>
rather than this sensor.
Read the
<a href=\"modelica://IBPSA.Fluid.Sensors.UsersGuide\">
IBPSA.Fluid.Sensors.UsersGuide</a>
prior to using this model to see about potential numerical problems if this sensor is used incorrectly
in a system model.
</p>
</html>",       revisions="<html>
<ul>
<li>
September 21, 2020, by Michael Wetter:<br/>
Introduced parameter <code>warnAboutOnePortConnection</code> and updated documentation.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1399\">#1399</a>.
</li>
<li>
February 21, 2020, by Michael Wetter:<br/>
Changed icon to display its operating state.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1294\">#1294</a>.
</li>
<li>
January 26, 2016 by Michael Wetter:<br/>
Added <code>quantity</code> attribute for mass fraction variables.<br/>
Made unit assignment of output signal final.
</li>
<li>
May 12, 2010 by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
      end RelativeHumidity;

      model RelativeHumidityTwoPort "Ideal two port relative humidity sensor"
        extends IBPSA.Fluid.Sensors.BaseClasses.PartialDynamicFlowSensor;
        extends Modelica.Icons.RoundSensor;
        Modelica.Blocks.Interfaces.RealOutput phi(final unit="1",
                                                  min=0,
                                                  start=phi_start)
          "Relative humidity of the passing fluid"
          annotation (Placement(transformation(extent={{-10,-10},{10,10}},
                rotation=90,
              origin={1,110})));
        parameter Real phi_start(final unit="1", min=0, max=1)=0.5
          "Initial or guess value of output (= state)"
          annotation (Dialog(group="Initialization"));

    protected
        Real phiMed(final unit="1", min=0, start=phi_start)
          "Relative humidity to which the sensor is exposed to";

    protected
        Modelica.Units.SI.Temperature T_a
          "Temperature of the medium flowing from port_a to port_b";
        Medium.MassFraction Xi_a[Medium.nXi](
          quantity=Medium.substanceNames[1:Medium.nXi])
          "Mass fraction of the medium flowing from port_a to port_b";
        Real phi_a(final unit="1")
          "Relative humidity of the medium flowing from port_a to port_b";
        Modelica.Units.SI.Temperature T_b
          "Temperature of the medium flowing from port_b to port_a";
        Medium.MassFraction Xi_b[Medium.nXi](
          quantity=Medium.substanceNames[1:Medium.nXi])
          "Mass fraction of the medium flowing from port_b to port_a";
        Real phi_b(final unit="1")
          "Relative humidity of the medium flowing from port_b to port_a";

      initial equation
        if dynamic then
          if initType == Modelica.Blocks.Types.Init.SteadyState then
            der(phi) = 0;
          elseif initType == Modelica.Blocks.Types.Init.InitialState or
                 initType == Modelica.Blocks.Types.Init.InitialOutput then
            phi = phi_start;
          end if;
        end if;
      equation

        // Since the sensor does not affect the medium, we can use
        // port_b.Xi_outflow = inStream(port_a.Xi_outflow).
        Xi_a = port_b.Xi_outflow;

        T_a=Medium.temperature_phX(
            p=port_a.p,
            h=port_b.h_outflow,
            X=Xi_a);

        phi_a = IBPSA.Utilities.Psychrometrics.Functions.phi_pTX(
          p=port_a.p,
          T=T_a,
          X_w=Xi_a[1]);

        if allowFlowReversal then
          phi_b = IBPSA.Utilities.Psychrometrics.Functions.phi_pTX(
            p=port_b.p,
            T=T_b,
            X_w=Xi_b[1]);
          T_b=Medium.temperature_phX(
            p=port_b.p,
            h=port_a.h_outflow,
            X=Xi_b);
          Xi_b = port_a.Xi_outflow;
          phiMed = Modelica.Fluid.Utilities.regStep(
                     x=port_a.m_flow,
                     y1=phi_a,
                     y2=phi_b,
                     x_small=m_flow_small);
        else
          phi_b = 0;
          T_b = 273.15;
          Xi_b = zeros(Medium.nXi);
          phiMed = phi_a;
        end if;
        // Output signal of sensor
        if dynamic then
          der(phi) = (phiMed-phi)*k*tauInv;
        else
          phi = phiMed;
        end if;
      annotation (defaultComponentName="senRelHum",
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={
              Text(
                extent={{102,124},{6,95}},
                textColor={0,0,0},
                textString="phi"),
              Line(points={{0,100},{0,70}}, color={0,0,127}),
              Line(points={{-100,0},{-70,0}}, color={0,128,255}),
              Line(points={{70,0},{100,0}}, color={0,128,255}),
              Text(
                extent={{-20,120},{-140,70}},
                textColor={0,0,0},
                textString=DynamicSelect("", String(phi, leftJustified=false, significantDigits=2)))}),
        Documentation(info="<html>
<p>
This model outputs the relative humidity of the fluid flowing from
<code>port_a</code> to <code>port_b</code>.
The sensor is ideal, i.e., it does not influence the fluid.
</p>
<p>
Note that this sensor can only be used with media that contain the variable <code>phi</code>,
which is typically the case for moist air models.
</p>
<p>
If the parameter <code>tau</code> is non-zero, then its output
is computed using a first order differential equation.
Setting <code>tau=0</code> is <i>not</i> recommend. See
<a href=\"modelica://IBPSA.Fluid.Sensors.UsersGuide\">
IBPSA.Fluid.Sensors.UsersGuide</a> for an explanation.
</p>
</html>",       revisions="<html>
<ul>
<li>
February 21, 2020, by Michael Wetter:<br/>
Changed icon to display its operating state.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1294\">#1294</a>.
</li>
<li>
January 26, 2016 by Michael Wetter:<br/>
Added <code>quantity</code> attribute for mass fraction variables.<br/>
Made unit assignment of output signal final.
</li>
<li>
January 18, 2016 by Filip Jorissen:<br/>
Using parameter <code>tauInv</code>
since this now exists in
<a href=\"modelica://IBPSA.Fluid.Sensors.BaseClasses.PartialDynamicFlowSensor\">IBPSA.Fluid.Sensors.BaseClasses.PartialDynamicFlowSensor</a>.
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/372\">#372</a>.
</li>
<li>
June 3, 2011 by Michael Wetter:<br/>
Revised implementation to add dynamics in such a way that
the time constant increases as the mass flow rate tends to zero.
This significantly improves the numerics.
</li>
<li>
Feb. 5, 2011 by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
      end RelativeHumidityTwoPort;

      package BaseClasses "Package with base classes for IBPSA.Fluid.Sensors"
        extends Modelica.Icons.BasesPackage;

        partial model PartialAbsoluteSensor
          "Partial component to model a sensor that measures a potential variable"

          replaceable package Medium=Modelica.Media.Interfaces.PartialMedium
            "Medium in the sensor"
            annotation (choices(
                choice(redeclare package Medium = IBPSA.Media.Air "Moist air"),
                choice(redeclare package Medium = IBPSA.Media.Water "Water"),
                choice(redeclare package Medium =
                    IBPSA.Media.Antifreeze.PropyleneGlycolWater (
                      property_T=293.15,
                      X_a=0.40)
                      "Propylene glycol water, 40% mass fraction")));

          parameter Boolean warnAboutOnePortConnection = true
           "Set to false to suppress warning about potential numerical issues, see IBPSA.Fluid.Sensors.UsersGuide for more information"
           annotation(HideResult=true);
          Modelica.Fluid.Interfaces.FluidPort_a port(redeclare package Medium=Medium, m_flow(min=0))
            annotation (Placement(transformation(
                origin={0,-100},
                extent={{-10,-10},{10,10}},
                rotation=90)));

      protected
          parameter String instanceName = getInstanceName() "Name of the instance";
        initial equation
          assert(not warnAboutOnePortConnection,
          "Sensor " + instanceName + " can lead to numerical problems if connected to a scalar fluid port.
  Only connect it to a vectorized fluid port, such as used in 'IBPSA.Fluid.MixingVolumes`.
  See IBPSA.Fluid.Sensors.UsersGuide for more information.
  To disable this warning, set 'warnAboutOnePortConnection = false' in "         + instanceName + ".",
          level=AssertionLevel.warning);
        equation
          port.m_flow = 0;
          port.h_outflow = 0;
          port.Xi_outflow = zeros(Medium.nXi);
          port.C_outflow = zeros(Medium.nC);
          annotation (Documentation(info="<html>
<p>
Partial component to model an absolute sensor.
The component can be used for pressure sensor models.
Use for other properties such as temperature or density is discouraged, because the enthalpy at the connector can have different meanings, depending on the connection topology. For these properties, use
<a href=\"modelica://IBPSA.Fluid.Sensors.BaseClasses.PartialFlowSensor\">
IBPSA.Fluid.Sensors.BaseClasses.PartialFlowSensor</a>.
</p>
</html>",
        revisions="<html>
<ul>
<li>
September 20, 2020, by Michael Wetter:<br/>
Introduced parameter <code>warnAboutOnePortConnection</code> and added assertion with level set to warning.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1399\"> #1399</a>.
</li>
<li>
January 18, 2019, by Jianjun Hu:<br/>
Limited the media choice.
See <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1050\">#1050</a>.
</li>
<li>
September 7, 2018, by Michael Wetter:<br/>
Changed
<code>port(redeclare package Medium=Medium, m_flow(min=0, max=0))</code>
to
<code>port(redeclare package Medium=Medium, m_flow(min=0))</code>
to avoid in Dymola 2019FD01 beta1 the message
\"port.m_flow has the range [0,0] - which is suspicious since the max-value should be above the min-value\"
which causes an error in pedantic mode.
Note that the MSL also uses only a <code>min</code> value.
</li>
<li>
March 22, 2017, by Filip Jorissen:<br/>
Set <code>m_flow(max=0)</code>.
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/687\">#687</a>.
</li>
<li>
February 12, 2011, by Michael Wetter:<br/>
First implementation.
Implementation is based on <code>Modelica.Fluid</code>.
</li>
</ul>
</html>"),
        Icon(
          graphics={
            Bitmap(
              visible = warnAboutOnePortConnection,
              extent={{-96,-82},{-64,-50}},
              fileName="modelica://IBPSA/Resources/Images/Fluid/Sensors/warningIcon.png")}));
        end PartialAbsoluteSensor;

        partial model PartialDynamicFlowSensor
          "Partial component to model sensors that measure flow properties using a dynamic model"
          extends PartialFlowSensor;

          parameter Modelica.Units.SI.Time tau(min=0) = 1
            "Time constant at nominal flow rate (use tau=0 for steady-state sensor, but see user guide for potential problems)";
          parameter Modelica.Blocks.Types.Init initType = Modelica.Blocks.Types.Init.InitialState
            "Type of initialization (InitialState and InitialOutput are identical)"
          annotation(Evaluate=true, Dialog(group="Initialization"));
      protected
          Real k(start=1)
            "Gain to take flow rate into account for sensor time constant";
          final parameter Boolean dynamic = tau > 1E-10 or tau < -1E-10
            "Flag, true if the sensor is a dynamic sensor"
            annotation(Evaluate=true);
          Real mNor_flow "Normalized mass flow rate";
          final parameter Real tauInv(final unit="s-1")= if dynamic then 1/tau else 0
            "Inverse of tau";
        equation
          if dynamic then
            mNor_flow = port_a.m_flow/m_flow_nominal;
            k = Modelica.Fluid.Utilities.regStep(x=port_a.m_flow,
                                                 y1= mNor_flow,
                                                 y2=-mNor_flow,
                                                 x_small=m_flow_small);
          else
            mNor_flow = 1;
            k = 1;
          end if;
          annotation (Icon(graphics={
                Line(visible=(tau > 0),
                points={{52,60},{58,74},{66,86},{76,92},{88,96},{98,96}}, color={0,
                      0,127})}), Documentation(info="<html>
<p>
Partial component to model a sensor that measures any intensive properties
of a flow, e.g., to get temperature or density in the flow
between fluid connectors.</p>
<p>
The sensor computes a gain that is zero at zero mass flow rate.
This avoids fast transients if the flow is close to zero, thereby
improving the numerical efficiency.
</p>
</html>",         revisions="<html>
<ul>
<li>
December 10, 2022, by Michael Wetter:<br/>
Corrected annotation to avoid comparing a real-valued parameter
for equality.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1671\">IBPSA, #1671</a>.
</li>
<li>
August 9, 2016, by Michael Wetter:<br/>
Improved documentation for <code>tau</code>.
</li>
<li>
January 12, 2016, by Filip Jorissen:<br/>
Added optional parameter <code>tauInv</code>.
</li>
<li>
May 29, 2014, by Michael Wetter:<br/>
Removed undesirable annotation <code>Evaluate=true</code>.
</li>
<li>
March 29, 2013, by Michael Wetter:<br/>
Changed the parameter <code>initType</code> to
<code>Modelica.Blocks.Types.Init.InitialState</code>.
This allows a pedantic model check in Dymola 2014 of models that instanciate sensors
but do not set this parameter. It also ensures that different Modelica simulators solve
the same initialization problem.
</li>
<li>
July 7, 2011, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
        end PartialDynamicFlowSensor;

        partial model PartialFlowSensor
          "Partial component to model sensors that measure flow properties"
          extends IBPSA.Fluid.Interfaces.PartialTwoPort;
          parameter Modelica.Units.SI.MassFlowRate m_flow_nominal(min=0)
            "Nominal mass flow rate, used for regularization near zero flow"
            annotation (Dialog(group="Nominal condition"));
          parameter Modelica.Units.SI.MassFlowRate m_flow_small(min=0) = 1E-4*
            m_flow_nominal
            "For bi-directional flow, temperature is regularized in the region |m_flow| < m_flow_small (m_flow_small > 0 required)"
            annotation (Dialog(tab="Advanced"));
        equation
          // mass balance
          port_b.m_flow = -port_a.m_flow;
          // momentum equation (no pressure loss)
          port_a.p = port_b.p;
          // isenthalpic state transformation (no storage and no loss of energy)
          port_a.h_outflow = if allowFlowReversal then inStream(port_b.h_outflow) else Medium.h_default;
          port_b.h_outflow = inStream(port_a.h_outflow);
          port_a.Xi_outflow = if allowFlowReversal then inStream(port_b.Xi_outflow) else Medium.X_default[1:Medium.nXi];
          port_b.Xi_outflow = inStream(port_a.Xi_outflow);
          port_a.C_outflow = if allowFlowReversal then inStream(port_b.C_outflow) else zeros(Medium.nC);
          port_b.C_outflow = inStream(port_a.C_outflow);
          annotation (Documentation(info="<html>
<p>
Partial component to model a sensor.
The sensor is ideal. It does not influence mass, energy,
species or substance balance, and it has no flow friction.
</p>
</html>",
        revisions="<html>
<ul>
<li>
August 15, 2015, by Filip Jorissen:<br/>
Implemented more efficient computation of <code>port_a.Xi_outflow</code>,
<code>port_a.h_outflow</code>
and <code>port_a.C_outflow</code> when <code>allowFlowReversal=false</code>.
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/281\">#281</a>.
</li>
<li>
June 19, 2015, by Michael Wetter:<br/>
Moved <code>m_flow_small</code> to the <code>Advanced</code> tab
as it usually need not be changed by the user.
Other models such as heat exchangers also have this parameter
on the <code>Advanced</code> tab.
</li>
<li>
February 12, 2011, by Michael Wetter:<br/>
First implementation.
Implementation is based on <code>Modelica.Fluid</code>.
</li>
</ul>
</html>"));
        end PartialFlowSensor;
      annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains base classes that are used to construct the models in
<a href=\"modelica://IBPSA.Fluid.Sensors\">IBPSA.Fluid.Sensors</a>.
</p>
</html>"));
      end BaseClasses;
    annotation (preferredView="info",
    Documentation(info="<html>
<p>
Package <code>Sensors</code> consists of idealized sensor components that
provide variables of a medium as
output signals. These signals can be, e.g., further processed
with components of the
<a href=\"modelica://Modelica.Blocks\">
Modelica.Blocks</a>
library.
</p>
</html>", revisions="<html>
<ul>
<li><i>22 Dec 2008</i>
    by R&uuml;diger Franke
    <ul>
    <li>flow sensors based on Modelica.Fluid.Interfaces.PartialTwoPort</li>
    <li>adapted documentation to stream connectors, i.e. less need for two port sensors</li>
    </ul>
</li>
<li><i>4 Dec 2008</i>
    by Michael Wetter<br/>
       included sensors for trace substance</li>
<li><i>31 Oct 2007</i>
    by Carsten Heinrich<br/>
       updated sensor models, included one and two port sensors for thermodynamic state variables</li>
</ul>
</html>"));
    end Sensors;

    package Interfaces "Package with interfaces for fluid models"
      extends Modelica.Icons.InterfacesPackage;

      partial model PartialTwoPort "Partial component with two ports"
        replaceable package Medium =
          Modelica.Media.Interfaces.PartialMedium "Medium in the component"
            annotation (choices(
              choice(redeclare package Medium = IBPSA.Media.Air "Moist air"),
              choice(redeclare package Medium = IBPSA.Media.Water "Water"),
              choice(redeclare package Medium =
                  IBPSA.Media.Antifreeze.PropyleneGlycolWater (
                    property_T=293.15,
                    X_a=0.40)
                    "Propylene glycol water, 40% mass fraction")));

        parameter Boolean allowFlowReversal = true
          "= false to simplify equations, assuming, but not enforcing, no flow reversal"
          annotation(Dialog(tab="Assumptions"), Evaluate=true);

        Modelica.Fluid.Interfaces.FluidPort_a port_a(
          redeclare final package Medium = Medium,
           m_flow(min=if allowFlowReversal then -Modelica.Constants.inf else 0),
           h_outflow(start = Medium.h_default, nominal = Medium.h_default))
          "Fluid connector a (positive design flow direction is from port_a to port_b)"
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b port_b(
          redeclare final package Medium = Medium,
          m_flow(max=if allowFlowReversal then +Modelica.Constants.inf else 0),
           h_outflow(start = Medium.h_default, nominal = Medium.h_default))
          "Fluid connector b (positive design flow direction is from port_a to port_b)"
          annotation (Placement(transformation(extent={{110,-10},{90,10}})));

        annotation (
          Documentation(info="<html>
<p>
This partial model defines an interface for components with two ports.
The treatment of the design flow direction and of flow reversal are predefined based on the parameter <code>allowFlowReversal</code>.
The component may transport fluid and may have internal storage for a given fluid <code>Medium</code>.
</p>
<h4>Implementation</h4>
<p>
This model is similar to
<a href=\"modelica://Modelica.Fluid.Interfaces.PartialTwoPort\">
Modelica.Fluid.Interfaces.PartialTwoPort</a>
but it does not use the <code>outer system</code> declaration.
This declaration is omitted as in building energy simulation,
many models use multiple media, an in practice,
users have not used this global definition to assign parameters.
</p>
</html>",       revisions="<html>
<ul>
<li>
January 18, 2019, by Jianjun Hu:<br/>
Limited the media choice.
See <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1050\">#1050</a>.
</li>
<li>
July 8, 2018, by Filip Jorissen:<br/>
Added nominal value of <code>h_outflow</code> in <code>FluidPorts</code>.
See <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/977\">#977</a>.
</li>
<li>
November 19, 2015, by Michael Wetter:<br/>
Removed parameters
<code>port_a_exposesState</code> and
<code>port_b_exposesState</code>
for <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/351\">#351</a>
and
<code>showDesignFlowDirection</code>
for <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/349\">#349</a>.
</li>
<li>
November 13, 2015, by Michael Wetter:<br/>
Assinged <code>start</code> attribute for leaving
enthalpy at <code>port_a</code> and <code>port_b</code>.
This was done to make the model similar to
<a href=\"modelica://IBPSA.Fluid.Interfaces.PartialFourPort\">
IBPSA.Fluid.Interfaces.PartialFourPort</a>.
</li>
<li>
November 12, 2015, by Michael Wetter:<br/>
Removed import statement.
</li>
<li>
October 21, 2014, by Michael Wetter:<br/>
Revised implementation.
Declared medium in ports to be <code>final</code>.
</li>
<li>
October 20, 2014, by Filip Jorisson:<br/>
First implementation.
</li>
</ul>
</html>"),Icon(coordinateSystem(
                preserveAspectRatio=true,
                extent={{-100,-100},{100,100}}), graphics={
              Polygon(
                points={{20,-70},{60,-85},{20,-100},{20,-70}},
                lineColor={0,128,255},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid,
                visible=not allowFlowReversal),
              Line(
                points={{55,-85},{-60,-85}},
                color={0,128,255},
                visible=not allowFlowReversal),
              Text(
                extent={{-149,-114},{151,-154}},
                textColor={0,0,255},
                textString="%name")}));
      end PartialTwoPort;
    annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains basic classes that are used to build
component models that change the state of the
fluid. The classes are not directly usable, but can
be extended when building a new model.
</p>
</html>"));
    end Interfaces;
  annotation (
  preferredView="info", Documentation(info="<html>
This package contains components for fluid flow systems such as
pumps, valves and sensors. For other fluid flow models, see
<a href=\"modelica://Modelica.Fluid\">Modelica.Fluid</a>.
</html>"),
  Icon(graphics={
          Polygon(points={{-70,26},{68,-44},{68,26},{2,-10},{-70,-42},{-70,26}},
              lineColor={0,0,0}),
          Line(points={{2,42},{2,-10}}),
          Rectangle(
            extent={{-18,50},{22,42}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid)}));
  end Fluid;

  package Media "Package with medium models"
    extends Modelica.Icons.Package;

    package Water "Package with model for liquid water with constant density"
       extends Modelica.Media.Water.ConstantPropertyLiquidWater(
         p_default=300000,
         reference_p=300000,
         reference_T=273.15,
         reference_X={1},
         AbsolutePressure(start=p_default),
         Temperature(start=T_default),
         Density(start=d_const),
         final cv_const=cp_const);
      // cp_const and cv_const have been made final because the model sets u=h.
      extends Modelica.Icons.Package;

      redeclare replaceable model BaseProperties "Base properties (p, d, T, h, u, R, MM and X and Xi) of a medium"
        parameter Boolean preferredMediumStates=false
          "= true if StateSelect.prefer shall be used for the independent property variables of the medium"
          annotation (Evaluate=true, Dialog(tab="Advanced"));
        final parameter Boolean standardOrderComponents=true
          "If true, and reducedX = true, the last element of X will be computed from the other ones";
        Modelica.Units.SI.Density d=d_const "Density of medium";
        Temperature T(stateSelect=
          if preferredMediumStates then StateSelect.prefer else StateSelect.default)
          "Temperature of medium";
        InputAbsolutePressure p "Absolute pressure of medium";
        InputMassFraction[nXi] Xi=fill(0, 0)
          "Structurally independent mass fractions";
        InputSpecificEnthalpy h "Specific enthalpy of medium";
        Modelica.Units.SI.SpecificInternalEnergy u
          "Specific internal energy of medium";

        Modelica.Units.SI.MassFraction[nX] X={1}
          "Mass fractions (= (component mass)/total mass  m_i/m)";
        final Modelica.Units.SI.SpecificHeatCapacity R_s=0
          "Gas constant (of mixture if applicable)";
        final Modelica.Units.SI.MolarMass MM=MM_const
          "Molar mass (of mixture or single fluid)";
        ThermodynamicState state
          "Thermodynamic state record for optional functions";


        Modelica.Units.NonSI.Temperature_degC T_degC=
            Modelica.Units.Conversions.to_degC(T) "Temperature of medium in [degC]";
        Modelica.Units.NonSI.Pressure_bar p_bar=Modelica.Units.Conversions.to_bar(p)
          "Absolute pressure of medium in [bar]";

        // Local connector definition, used for equation balancing check
        connector InputAbsolutePressure = input Modelica.Units.SI.AbsolutePressure
          "Pressure as input signal connector";
        connector InputSpecificEnthalpy = input Modelica.Units.SI.SpecificEnthalpy
          "Specific enthalpy as input signal connector";
        connector InputMassFraction = input Modelica.Units.SI.MassFraction
          "Mass fraction as input signal connector";

      equation
        h = cp_const*(T-reference_T);
        u = h;
        state.T = T;
        state.p = p;

        // Assertions to test for bounds
        assert(noEvent(T >= T_min), "In " + getInstanceName() + ": Temperature T = " + String(T) + " K exceeded its minimum allowed value of " +
      String(T_min-273.15) + " degC (" + String(T_min) + " Kelvin) as required from medium model \"IBPSA.Media.Water\".");

        assert(noEvent(T <= T_max), "In " + getInstanceName() + ": Temperature T = " + String(T) + " K exceeded its maximum allowed value of " +
      String(T_max-273.15) + " degC (" + String(T_max) + " Kelvin) as required from medium model \"IBPSA.Media.Water\".");

        assert(noEvent(p >= 0.0), "Pressure (= " + String(p) + " Pa) of medium \"IBPSA.Media.Water\" is negative\n(Temperature = " + String(T) + " K)");

        annotation(Documentation(info="<html>
<p>
Model with basic thermodynamic properties.
</p>
<p>
This base properties model is identical to
<a href=\"modelica://Modelica.Media.Water.ConstantPropertyLiquidWater\">
Modelica.Media.Water.ConstantPropertyLiquidWater</a>,
except that the equation
<code>u = cv_const*(T - reference_T)</code>
has been replaced by <code>u=h</code> because
<code>cp_const=cv_const</code>.
</p>
<p>
This model provides equation for the following thermodynamic properties:
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\" summary=\"Thermodynamic properties\">
  <tr><td><strong>Variable</strong></td>
      <td><strong>Unit</strong></td>
      <td><strong>Description</strong></td></tr>
  <tr><td>T</td>
      <td>K</td>
      <td>temperature</td></tr>
  <tr><td>p</td>
      <td>Pa</td>
      <td>absolute pressure</td></tr>
  <tr><td>d</td>
      <td>kg/m3</td>
      <td>density</td></tr>
  <tr><td>h</td>
      <td>J/kg</td>
      <td>specific enthalpy</td></tr>
  <tr><td>u</td>
      <td>J/kg</td>
      <td>specific internal energy</td></tr>
  <tr><td>Xi[nXi]</td>
      <td>kg/kg</td>
      <td>independent mass fractions m_i/m</td></tr>
  <tr><td>R</td>
      <td>J/kg.K</td>
      <td>gas constant</td></tr>
  <tr><td>M</td>
      <td>kg/mol</td>
      <td>molar mass</td></tr>
</table>
</html>"));
      end BaseProperties;

    function enthalpyOfLiquid "Return the specific enthalpy of liquid"
      extends Modelica.Icons.Function;
        input Modelica.Units.SI.Temperature T "Temperature";
        output Modelica.Units.SI.SpecificEnthalpy h "Specific enthalpy";
    algorithm
      h := cp_const*(T-reference_T);
    annotation (
      smoothOrder=5,
      Inline=true,
    Documentation(info="<html>
<p>
Enthalpy of the water.
</p>
</html>",     revisions="<html>
<ul>
<li>
October 16, 2014 by Michael Wetter:<br/>
First implementation.
This function is used by
<a href=\"modelica://IBPSA.Fluid.MixingVolumes.MixingVolumeMoistAir\">
IBPSA.Fluid.MixingVolumes.MixingVolumeMoistAir</a>.
</li>
</ul>
</html>"));
    end enthalpyOfLiquid;
      annotation(Documentation(info="<html>
<p>
This medium package models liquid water.
</p>
<p>
The mass density is computed using a constant value of <i>995.586</i> kg/s.
For a medium model in which the density is a function of temperature, use
<a href=\"modelica://IBPSA.Media.Specialized.Water.TemperatureDependentDensity\">
IBPSA.Media.Specialized.Water.TemperatureDependentDensity</a> which may have considerably higher computing time.
</p>
<p>
For the specific heat capacities at constant pressure and at constant volume,
a constant value of <i>4184</i> J/(kg K), which corresponds to <i>20</i>&deg;C
is used.
The figure below shows the relative error of the specific heat capacity that
is introduced by this simplification.
</p>
<p align=\"center\">
<img src=\"modelica://IBPSA/Resources/Images/Media/Water/plotCp.png\" border=\"1\"
alt=\"Relative variation of specific heat capacity with temperature\"/>
</p>
<p>
The enthalpy is computed using the convention that <i>h=0</i>
if <i>T=0</i> &deg;C.
</p>
<h4>Limitations</h4>
<p>
Density, specific heat capacity, thermal conductivity and viscosity are constant.
Water is modeled as an incompressible liquid.
There are no phase changes.
</p>
</html>",     revisions="<html>
<ul>
<li>
September 28, 2020, by Michael Wetter:<br/>
Reformulated <code>BaseProperties</code> to avoid event-triggering assertions.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1401\">#1401</a>.
</li>
<li>
October 26, 2018, by Filip Jorissen and Michael Wetter:<br/>
Now printing different messages if temperature is above or below its limit,
and adding instance name as JModelica does not print the full instance name in the assertion.
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1045\">#1045</a>.
</li>
<li>
June 6, 2015, by Michael Wetter:<br/>
Set <code>AbsolutePressure(start=p_default)</code> to avoid
a translation error if
<a href=\"modelica://IBPSA.Fluid.Sources.Examples.TraceSubstancesFlowSource\">
IBPSA.Fluid.Sources.Examples.TraceSubstancesFlowSource</a>
(if used with water instead of air)
is translated in pedantic mode in Dymola 2016.
The reason is that pressures use <code>Medium.p_default</code> as start values,
but
<a href=\"modelica://Modelica.Media.Interfaces.Types\">
Modelica.Media.Interfaces.Types</a>
sets a default value of <i>1E-5</i>.
A similar change has been done for pressure and density.
This fixes
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/266\">#266</a>.
</li>
<li>
June 6, 2015, by Michael Wetter:<br/>
Changed type of <code>BaseProperties.T</code> from
<code>Modelica.Units.SI.Temperature</code> to <code>Temperature</code>.
Otherwise, it has a different start value than <code>Medium.T</code>, which
causes an error if
<a href=\"modelica://IBPSA.Media.Examples.WaterProperties\">
IBPSA.Media.Examples.WaterProperties</a>
is translated in pedantic mode.
This fixes
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/266\">#266</a>.
</li>
<li>
June 5, 2015, by Michael Wetter:<br/>
Added <code>stateSelect</code> attribute in <code>BaseProperties.T</code>
to allow correct use of <code>preferredMediumState</code> as
described in
<a href=\"modelica://Modelica.Media.Interfaces.PartialMedium\">
Modelica.Media.Interfaces.PartialMedium</a>,
and set <code>preferredMediumState=false</code>
to keep the same states as were used before.
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/260\">#260</a>.
</li>
<li>
June 5, 2015, by Michael Wetter:<br/>
Removed <code>ThermodynamicState</code> declaration as this lead to
the error
\"Attempting to redeclare record ThermodynamicState when the original was not replaceable.\"
in Dymola 2016 using the pedantic model check.
</li>
<li>
May 1, 2015, by Michael Wetter:<br/>
Added <code>Inline=true</code> for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/227\">
issue 227</a>.
</li>
<li>
February 25, 2015, by Michael Wetter:<br/>
Removed <code>stateSelect</code> attribute on pressure as this caused
<a href=\"modelica://IBPSA.Examples.Tutorial.SpaceCooling.System3\">
IBPSA.Examples.Tutorial.SpaceCooling.System3</a>
to fail with the error message
\"differentiated if-then-else was not continuous\".
</li>
<li>
October 15, 2014, by Michael Wetter:<br/>
Reimplemented media based on
<a href=\"https://github.com/ibpsa/modelica-ibpsa/blob/446aa83720884052476ad6d6d4f90a6a29bb8ec9/IBPSA/Media/Water.mo\">446aa83</a>.
</li>
<li>
November 15, 2013, by Michael Wetter:<br/>
Complete new reimplementation because the previous version
had the option to add a compressibility to the medium, which
has never been used.
</li>
</ul>
</html>"),
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
            graphics={
            Polygon(
              points={{16,-28},{32,-42},{26,-48},{10,-36},{16,-28}},
              lineColor={95,95,95},
              fillPattern=FillPattern.Sphere,
              fillColor={95,95,95}),
            Polygon(
              points={{10,34},{26,44},{30,36},{14,26},{10,34}},
              lineColor={95,95,95},
              fillPattern=FillPattern.Sphere,
              fillColor={95,95,95}),
            Ellipse(
              extent={{-82,52},{24,-54}},
              lineColor={95,95,95},
              fillPattern=FillPattern.Sphere,
              fillColor={0,0,0}),
            Ellipse(
              extent={{22,82},{80,24}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Sphere,
              fillColor={95,95,95}),
            Ellipse(
              extent={{20,-30},{78,-88}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Sphere,
              fillColor={95,95,95})}));
    end Water;
    annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains media models for water and moist air.
The media models in this package are
compatible with
<a href=\"modelica://Modelica.Media\">
Modelica.Media</a>
but the implementation is in general simpler, which often
leads to more efficient simulation.
Due to the simplifications, the media model of this package
are generally accurate for a smaller temperature range than the
models in <a href=\"modelica://Modelica.Media\">
Modelica.Media</a>, but the smaller temperature range may often be
sufficient for building HVAC applications.
</p>
</html>"),
  Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
          graphics={
          Line(
            points = {{-76,-80},{-62,-30},{-32,40},{4,66},{48,66},{73,45},{62,-8},{48,-50},{38,-80}},
            color={64,64,64},
            smooth=Smooth.Bezier),
          Line(
            points={{-40,20},{68,20}},
            color={175,175,175}),
          Line(
            points={{-40,20},{-44,88},{-44,88}},
            color={175,175,175}),
          Line(
            points={{68,20},{86,-58}},
            color={175,175,175}),
          Line(
            points={{-60,-28},{56,-28}},
            color={175,175,175}),
          Line(
            points={{-60,-28},{-74,84},{-74,84}},
            color={175,175,175}),
          Line(
            points={{56,-28},{70,-80}},
            color={175,175,175}),
          Line(
            points={{-76,-80},{38,-80}},
            color={175,175,175}),
          Line(
            points={{-76,-80},{-94,-16},{-94,-16}},
            color={175,175,175})}));
  end Media;

  package Utilities "Package with utility functions such as for I/O"
    extends Modelica.Icons.Package;

    package Math "Library with functions such as for smoothing"
      extends Modelica.Icons.Package;

      package Functions "Package with mathematical functions"
        extends Modelica.Icons.VariantsPackage;

        function regStep
          "Approximation of a general step, such that the approximation is continuous and differentiable"
          extends Modelica.Icons.Function;
          input Real x "Abscissa value";
          input Real y1 "Ordinate value for x > 0";
          input Real y2 "Ordinate value for x < 0";
          input Real x_small(min=0) = 1e-5
            "Approximation of step for -x_small <= x <= x_small; x_small >= 0 required";
          output Real y "Ordinate value to approximate y = if x > 0 then y1 else y2";
        algorithm
          y := smooth(1, if x >  x_small then y1 else
                         if x < -x_small then y2 else
                         if x_small > 0 then (x/x_small)*((x/x_small)^2 - 3)*(y2-y1)/4 + (y1+y2)/2 else (y1+y2)/2);

          annotation(Inline=true,
          Documentation(revisions="<html>
<ul>
<li><i>February 18, 2016</i>
    by Marcus Fuchs:<br/>
    Add function with <code>Inline = true</code> in annotations to package for better performance,
    as suggested in <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/300\">#300</a> .</li>
<li><i>August 12, 2008</i>
    by <a href=\"mailto:Michael.Sielemann@dlr.de\">Michael Sielemann</a>:<br/>
    Minor modification to cover the limit case <code>x_small -> 0</code> without division by zero.</li>
<li><i>April 29, 2008</i>
    by <a href=\"mailto:Martin.Otter@DLR.de\">Martin Otter</a>:<br/>
    Designed and implemented.</li>
</ul>
</html>",         info="<html>
<p>
This function is used to approximate the equation
</p>
<pre>
    y = <b>if</b> x &gt; 0 <b>then</b> y1 <b>else</b> y2;
</pre>

<p>
by a smooth characteristic, so that the expression is continuous and differentiable:
</p>

<pre>
   y = <b>smooth</b>(1, <b>if</b> x &gt;  x_small <b>then</b> y1 <b>else</b>
                 <b>if</b> x &lt; -x_small <b>then</b> y2 <b>else</b> f(y1, y2));
</pre>

<p>
In the region <code>-x_small &lt; x &lt; x_small</code> a 2nd order polynomial is used
for a smooth transition from <code>y1</code> to <code>y2</code>.
</p>
</html>"));
        end regStep;
      annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains functions for commonly used
mathematical operations. The functions are used in
the blocks
<a href=\"modelica://IBPSA.Utilities.Math\">
IBPSA.Utilities.Math</a>.
</p>
</html>"));
      end Functions;
    annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains blocks and functions for commonly used
mathematical operations.
The classes in this package augment the classes
<a href=\"modelica://Modelica.Blocks\">
Modelica.Blocks</a>.
</p>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
              {100,100}}), graphics={Line(points={{-80,0},{-68.7,34.2},{-61.5,53.1},
                {-55.1,66.4},{-49.4,74.6},{-43.8,79.1},{-38.2,79.8},{-32.6,76.6},{
                -26.9,69.7},{-21.3,59.4},{-14.9,44.1},{-6.83,21.2},{10.1,-30.8},{17.3,
                -50.2},{23.7,-64.2},{29.3,-73.1},{35,-78.4},{40.6,-80},{46.2,-77.6},
                {51.9,-71.5},{57.5,-61.9},{63.9,-47.2},{72,-24.8},{80,0}}, color={
                0,0,0}, smooth=Smooth.Bezier)}));
    end Math;

    package Psychrometrics "Library with psychrometric functions"
      extends Modelica.Icons.VariantsPackage;

      package Constants "Library of constants for psychometric functions"
        extends Modelica.Icons.Package;

        constant Real k_mair = 0.6219647130774989 "Ratio of molar weights";
        annotation (
          Documentation(info="<html>
<p>
This package provides constants for functions used
in the calculation of thermodynamic properties of moist air.
</p>
</html>",       revisions="<html>
<ul>
<li>
May 24, 2016, by Filip Jorissen:<br/>
Added reference temperature.
</li>
<li>
July 24, 2014, by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"),Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
            Polygon(
              origin={-9.2597,25.6673},
              fillColor={102,102,102},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{48.017,11.336},{48.017,11.336},{10.766,11.336},{-25.684,10.95},{-34.944,-15.111},{-34.944,-15.111},{-32.298,-15.244},{-32.298,-15.244},{-22.112,0.168},{11.292,0.234},{48.267,-0.097},{48.267,-0.097}},
              smooth=Smooth.Bezier),
            Polygon(
              origin={-19.9923,-8.3993},
              fillColor={102,102,102},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{3.239,37.343},{3.305,37.343},{-0.399,2.683},{-16.936,-20.071},{-7.808,-28.604},{6.811,-22.519},{9.986,37.145},{9.986,37.145}},
              smooth=Smooth.Bezier),
            Polygon(
              origin={23.753,-11.5422},
              fillColor={102,102,102},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-10.873,41.478},{-10.873,41.478},{-14.048,-4.162},{-9.352,-24.8},{7.912,-24.469},{16.247,0.27},{16.247,0.27},{13.336,0.071},{13.336,0.071},{7.515,-9.983},{-3.134,-7.271},{-2.671,41.214},{-2.671,41.214}},
              smooth=Smooth.Bezier)}));
      end Constants;

      package Functions "Package with psychrometric functions"
        extends Modelica.Icons.Package;

        function X_pTphi
          "Absolute humidity for given pressure, dry bulb temperature and relative humidity"
          extends Modelica.Icons.Function;
          input Modelica.Units.SI.Pressure p "Absolute pressure of the medium";
          input Modelica.Units.SI.Temperature T "Dry bulb temperature";
          input Real phi(unit="1") "Relative humidity";
          output Modelica.Units.SI.MassFraction X_w
            "Water vapor mass fraction per unit mass total air";

        algorithm
          X_w:=phi/((p/saturationPressure(T)-phi) / IBPSA.Utilities.Psychrometrics.Constants.k_mair + phi);
          annotation (
            inverse(phi=phi_pTX(p,T,X_w)),
            smoothOrder=1,
            Documentation(info="<html>
<p>
Absolute humidity of air for given
pressure, temperature and relative humidity.
</p>
<p>
Note that the absolute humidity is in <i>kg/kg</i>
total air, and not dry air.
</p>
</html>",
        revisions="<html>
<ul>
<li>
April 4, 2019 by Filip Jorissen:<br/>
First implementation.
</li>
</ul>
</html>"));
        end X_pTphi;

        function phi_pTX
          "Relative humidity for given pressure, dry bulb temperature and moisture mass fraction"
          extends Modelica.Icons.Function;
          input Modelica.Units.SI.Pressure p "Absolute pressure of the medium";
          input Modelica.Units.SI.Temperature T "Dry bulb temperature";
          input Modelica.Units.SI.MassFraction X_w
            "Water vapor mass fraction per unit mass total air";
          output Real phi(unit="1") "Relative humidity";
        algorithm
          phi :=p/saturationPressure(T)*X_w/(X_w +
            IBPSA.Utilities.Psychrometrics.Constants.k_mair*(1-X_w));
          annotation (
            inverse(X_w=X_pTphi(p,T,phi)),
            smoothOrder=1,
            Documentation(info="<html>
<p>
Relative humidity of air for given
pressure, temperature and water vapor mass fraction.
</p>
<p>
Note that the water vapor mass fraction must be in <i>kg/kg</i>
total air, and not dry air.
</p>
</html>",
        revisions="<html>
<ul>
<li>
April 4, 2019 by Filip Jorissen:<br/>
Added inverse annotation
for <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1110\">#1110</a>.
</li>
<li>
November 17, 2014 by Michael Wetter:<br/>
Removed test that constrains the saturation pressure to be
lower than <code>p</code>.
I do not see any numerical problems without this test.
</li>
<li>
November 13, 2014 by Michael Wetter:<br/>
First implementation.
</li>
</ul>
</html>"));
        end phi_pTX;

        function saturationPressure
          "Saturation curve valid for 223.16 <= T <= 373.16 (and slightly outside with less accuracy)"
          extends Modelica.Icons.Function;
          input Modelica.Units.SI.Temperature TSat(displayUnit="degC", nominal=300)
            "Saturation temperature";
          output Modelica.Units.SI.AbsolutePressure pSat(displayUnit="Pa", nominal=1000)
            "Saturation pressure";

        algorithm
          pSat := IBPSA.Utilities.Math.Functions.regStep(
                     y1=IBPSA.Utilities.Psychrometrics.Functions.saturationPressureLiquid(TSat),
                     y2=IBPSA.Utilities.Psychrometrics.Functions.sublimationPressureIce(TSat),
                     x=TSat-273.16,
                     x_small=1.0);
          annotation(Inline=true,
            smoothOrder=1,
            Documentation(info="<html>
<p>
Saturation pressure of water, computed from temperature,
according to Wagner <i>et al.</i> (1993).
The range of validity is between
<i>190</i> and <i>373.16</i> Kelvin.
</p>
<h4>References</h4>
<p>
Wagner W., A. Saul, A. Pruss.
 <i>International equations for the pressure along the melting and along the sublimation curve of ordinary water substance</i>,
equation 3.5. 1993.
<a href=\"http://aip.scitation.org/doi/pdf/10.1063/1.555947?class=pdf\">
http://aip.scitation.org/doi/pdf/10.1063/1.555947?class=pdf</a>.
</p>
</html>",
        revisions="<html>
<ul>
<li>
March 15, 2016, by Michael Wetter:<br/>
Replaced <code>spliceFunction</code> with <code>regStep</code>.
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/300\">issue 300</a>.
</li>
<li>
August 19, 2015 by Michael Wetter:<br/>
Changed <code>smoothOrder</code> from <i>5</i> to <i>1</i> as
<a href=\"modelica://IBPSA.Utilities.Math.Functions.spliceFunction\">
IBPSA.Utilities.Math.Functions.spliceFunction</a> is only once
continuously differentiable.
Inlined the function.
</li>
<li>
November 20, 2013 by Michael Wetter:<br/>
First implementation, moved from <code>IBPSA.Media</code>.
</li>
</ul>
</html>"));
        end saturationPressure;

        function saturationPressureLiquid
          "Return saturation pressure of water as a function of temperature T in the range of 273.16 to 373.16 K"
          extends Modelica.Icons.Function;
          input Modelica.Units.SI.Temperature TSat(displayUnit="degC", nominal=300)
            "Saturation temperature";
          output Modelica.Units.SI.AbsolutePressure pSat(displayUnit="Pa", nominal=1000)
            "Saturation pressure";
        algorithm
          pSat := 611.657*Modelica.Math.exp(17.2799 - 4102.99/(TSat - 35.719));

          annotation (
            smoothOrder=99,
            derivative=IBPSA.Utilities.Psychrometrics.Functions.BaseClasses.der_saturationPressureLiquid,
            Inline=true,
            Documentation(info="<html>
<p>
Saturation pressure of water above the triple point temperature computed from temperature
according to Wagner <i>et al.</i> (1993). The range of validity is between
<i>273.16</i> and <i>373.16</i> Kelvin.
</p>
<h4>References</h4>
<p>
Wagner W., A. Saul, A. Pruss.
 <i>International equations for the pressure along the melting and along the sublimation curve of ordinary water substance</i>,
equation 3.5. 1993.
<a href=\"http://aip.scitation.org/doi/pdf/10.1063/1.555947?class=pdf\">
http://aip.scitation.org/doi/pdf/10.1063/1.555947?class=pdf</a>.
</p>
</html>",
        revisions="<html>
<ul>
<li>
November 20, 2013 by Michael Wetter:<br/>
First implementation, moved from <code>IBPSA.Media</code>.
</li>
</ul>
</html>"));
        end saturationPressureLiquid;

        function sublimationPressureIce
          "Return sublimation pressure of water as a function of temperature T between 190 and 273.16 K"
          extends Modelica.Icons.Function;
          input Modelica.Units.SI.Temperature TSat(displayUnit="degC", nominal=300)
            "Saturation temperature";
          output Modelica.Units.SI.AbsolutePressure pSat(displayUnit="Pa", nominal=1000)
            "Saturation pressure";
      protected
          Modelica.Units.SI.Temperature TTriple=273.16 "Triple point temperature";
          Modelica.Units.SI.AbsolutePressure pTriple=611.657 "Triple point pressure";
          Real r1=TSat/TTriple "Common subexpression";
          Real a[2]={-13.9281690,34.7078238} "Coefficients a[:]";
          Real n[2]={-1.5,-1.25} "Coefficients n[:]";
        algorithm
          pSat := exp(a[1] - a[1]*r1^n[1] + a[2] - a[2]*r1^n[2])*pTriple;
          annotation (
            Inline=false,
            smoothOrder=5,
            derivative=IBPSA.Utilities.Psychrometrics.Functions.BaseClasses.der_sublimationPressureIce,
            Documentation(info="<html>
<p>
Sublimation pressure of water below the triple point temperature, computed from temperature,
according to Wagner <i>et al.</i> (1993).
The range of validity is between
<i>190</i> and <i>273.16</i> Kelvin.
</p>
<h4>References</h4>
<p>
Wagner W., A. Saul, A. Pruss.
 <i>International equations for the pressure along the melting and along the sublimation curve of ordinary water substance</i>,
equation 3.5. 1993.
<a href=\"http://aip.scitation.org/doi/pdf/10.1063/1.555947?class=pdf\">
http://aip.scitation.org/doi/pdf/10.1063/1.555947?class=pdf</a>.
</p>
</html>",
        revisions="<html>
<ul>
<li>
November 20, 2013 by Michael Wetter:<br/>
First implementation, moved from <code>IBPSA.Media</code>.
</li>
</ul>
</html>"));
        end sublimationPressureIce;

        package BaseClasses "Package with base classes for IBPSA.Utilities.Psychrometrics.Functions"
          extends Modelica.Icons.BasesPackage;

          function der_saturationPressureLiquid
            "Derivative of the function saturationPressureLiquid"
            extends Modelica.Icons.Function;
            input Modelica.Units.SI.Temperature TSat "Saturation temperature";
            input Real dTSat(unit="K/s") "Saturation temperature derivative";
            output Real psat_der(unit="Pa/s") "Differential of saturation pressure";

          algorithm
            psat_der:=611.657*Modelica.Math.exp(17.2799 - 4102.99
                      /(TSat - 35.719))*4102.99*dTSat/(TSat - 35.719)^2;

            annotation(Inline=false,
              smoothOrder=98,
              Documentation(info="<html>
<p>
Derivative of function
<a href=\"modelica://IBPSA.Utilities.Psychrometrics.Functions.saturationPressureLiquid\">
IBPSA.Utilities.Psychrometrics.Functions.saturationPressureLiquid</a>.
</p>
</html>", revisions="<html>
<ul>
<li>
September 12, 2020, by Michael Wetter:<br/>
Corrected name of argument to comply with derivative specification.<br/>
This is for
<a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1393\">#1393</a>.
</li>
<li>
November 20, 2013 by Michael Wetter:<br/>
First implementation, moved from <code>IBPSA.Media</code>.
</li>
</ul>
</html>"));
          end der_saturationPressureLiquid;

          function der_sublimationPressureIce
            "Derivative of function sublimationPressureIce"
              extends Modelica.Icons.Function;
            input Modelica.Units.SI.Temperature TSat(displayUnit="degC", nominal=300)
              "Saturation temperature";
              input Real dTSat(unit="K/s") "Sublimation temperature derivative";
              output Real psat_der(unit="Pa/s") "Sublimation pressure derivative";
        protected
            Modelica.Units.SI.Temperature TTriple=273.16 "Triple point temperature";
            Modelica.Units.SI.AbsolutePressure pTriple=611.657 "Triple point pressure";
              Real r1=TSat/TTriple "Common subexpression 1";
              Real r1_der=dTSat/TTriple "Derivative of common subexpression 1";
              Real a[2]={-13.9281690,34.7078238} "Coefficients a[:]";
              Real n[2]={-1.5,-1.25} "Coefficients n[:]";
          algorithm
              psat_der := exp(a[1] - a[1]*r1^n[1] + a[2] - a[2]*r1^n[2])*pTriple*(-(a[1]
                *(r1^(n[1] - 1)*n[1]*r1_der)) - (a[2]*(r1^(n[2] - 1)*n[2]*r1_der)));
              annotation (
                Inline=false,
                smoothOrder=5,
                Documentation(info="<html>
<p>
Derivative of function
<a href=\"modelica://IBPSA.Utilities.Psychrometrics.Functions.sublimationPressureIce\">
IBPSA.Utilities.Psychrometrics.Functions.sublimationPressureIce</a>.
</p>
</html>", revisions="<html>
<ul>
<li>
September 12, 2020, by Michael Wetter:<br/>
Change name of argument <code>dTsat</code> to <code>dTSat</code> for consistency
with
<a href=\"modelica://IBPSA.Utilities.Psychrometrics.Functions.BaseClasses.der_saturationPressureLiquid\">
IBPSA.Utilities.Psychrometrics.Functions.BaseClasses.der_saturationPressureLiquid</a>.
</li>
<li>
November 20, 2013 by Michael Wetter:<br/>
First implementation, moved from <code>IBPSA.Media</code>.
</li>
</ul>
</html>"));
          end der_sublimationPressureIce;
        annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains base classes that are used to construct the models in
<a href=\"modelica://IBPSA.Utilities.Psychrometrics.Functions\">IBPSA.Utilities.Psychrometrics.Functions</a>.
</p>
</html>"));
        end BaseClasses;
        annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains functions for psychrometric calculations.
</p>

The nomenclature used in this package is described at
<a href=\"modelica://IBPSA.UsersGuide.Conventions\">
IBPSA.UsersGuide.Conventions</a>.
</html>"));
      end Functions;
    annotation (preferredView="info", Documentation(info="<html>
<p>
This package contains blocks and functions for psychrometric calculations.
</p>
<p>
The nomenclature used in this package is described at
<a href=\"modelica://IBPSA.UsersGuide.Conventions\">
IBPSA.UsersGuide.Conventions</a>.
</p>
</html>"));
    end Psychrometrics;
  annotation (
  preferredView="info", Documentation(info="<html>
<p>
This package contains utility models such as for thermal comfort calculation, input/output, co-simulation, psychrometric calculations and various functions that are used throughout the library.
</p>
</html>"),
  Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
      Polygon(
        origin={1.3835,-4.1418},
        rotation=45.0,
        fillColor={64,64,64},
        pattern=LinePattern.None,
        fillPattern=FillPattern.Solid,
        points={{-15.0,93.333},{-15.0,68.333},{0.0,58.333},{15.0,68.333},{15.0,93.333},{20.0,93.333},{25.0,83.333},{25.0,58.333},{10.0,43.333},{10.0,-41.667},{25.0,-56.667},{25.0,-76.667},{10.0,-91.667},{0.0,-91.667},{0.0,-81.667},{5.0,-81.667},{15.0,-71.667},{15.0,-61.667},{5.0,-51.667},{-5.0,-51.667},{-15.0,-61.667},{-15.0,-71.667},{-5.0,-81.667},{0.0,-81.667},{0.0,-91.667},{-10.0,-91.667},{-25.0,-76.667},{-25.0,-56.667},{-10.0,-41.667},{-10.0,43.333},{-25.0,58.333},{-25.0,83.333},{-20.0,93.333}}),
      Polygon(
        origin={10.1018,5.218},
        rotation=-45.0,
        fillColor={255,255,255},
        fillPattern=FillPattern.Solid,
        points={{-15.0,87.273},{15.0,87.273},{20.0,82.273},{20.0,27.273},{10.0,17.273},{10.0,7.273},{20.0,2.273},{20.0,-2.727},{5.0,-2.727},{5.0,-77.727},{10.0,-87.727},{5.0,-112.727},{-5.0,-112.727},{-10.0,-87.727},{-5.0,-77.727},{-5.0,-2.727},{-20.0,-2.727},{-20.0,2.273},{-10.0,7.273},{-10.0,17.273},{-20.0,27.273},{-20.0,82.273}})}));
  end Utilities;
annotation (
version="4.0.0",
versionDate="2018-09-28",
dateModified = "2018-09-28",
uses(Modelica(version="4.0.0")),
conversion(from(version="3.0.0",
                script="modelica://IBPSA/Resources/Scripts/Conversion/ConvertIBPSA_from_3.0_to_4.0.mos")),
preferredView="info",
Documentation(info="<html>
<p>
<img
align=\"right\"
alt=\"Logo of IBPSA\"
src=\"modelica://IBPSA/Resources/Images/IBPSA-logo-text.png\" border=\"1\"/>
The <code>IBPSA</code> library is a free library
that provides more than 300 classes (models, functions, etc.) for the development of
Modelica libraries for building and community energy and control systems.
The library is compatible with models from the Modelica Standard Library,
in particular with models from
<a href=\"modelica://Modelica.Fluid\">Modelica.Fluid</a>
and
<a href=\"modelica://Modelica.Media\">Modelica.Media</a>.
</p>
<p>
The development of the IBPSA library is organized through the
<a href=\"https://ibpsa.github.io/project1\">IBPSA Project 1</a>
of the International Building Performance Simulation Association (IBPSA).
From 2012 to 2017, the development was organized through the
<a href=\"http://www.iea-annex60.org\">Annex 60 project</a>
of the Energy in Buildings and Communities Programme of the International Energy Agency (IEA EBC).
</p>
<p>
The intent of the library is that it will be extended by
implementations of Modelica libraries that are targeted to end-users.
Major goals are
</p>
<ul>
<li>to codify best practice and to provide a solid foundation onto which
other libraries for building and community energy systems can be built, and
</li>
<li>
to avoid a fragmentation of libraries that serve similar purpose but
that cannot share models among each others, thereby duplicating efforts
for model development and validation.
</li>
</ul>
<p>
Hence, this library is typically not used directly by end-users,
but rather by developers of libraries that will be distributed to end-users.
Libraries that are using the <code>IBPSA</code> library as their core, or
that are working on using the <code>IBPSA</code> as their core, include, in
alphabetic order:
</p>
<ul>
<li>
The <code>AixLib</code> library from RWTH Aachen, Germany, available at
<a href=\"https://github.com/RWTH-EBC/AixLib\">https://github.com/RWTH-EBC/AixLib</a>
</li>
<li>
The <code>Buildings</code> library from Lawrence Berkeley National Laboratory, Berkeley, CA, available at
<a href=\"http://simulationresearch.lbl.gov/modelica\">http://simulationresearch.lbl.gov/modelica/</a>.
</li>
<li>
The <code>BuildingSystems</code> library from
Universit&auml;t der K&uuml;nste Berlin, Germany,
available at
<a href=\"http://www.modelica-buildingsystems.de/\">http://www.modelica-buildingsystems.de/</a>.
</li>
<li>
The <code>IDEAS</code> library from KU Leuven, Belgium, available at
<a href=\"https://github.com/open-ideas/IDEAS\">https://github.com/open-ideas/IDEAS</a>.
</li>
</ul>
<p>
The library also contains more than 300 example and validation models. For Dymola,
each of these example and validation models contains a script that simulates it and
plots certain variables of interest.
</p>
<p>
The web page for this library is
<a href=\"https://github.com/ibpsa/modelica\">https://github.com/ibpsa/modelica</a>.
Contributions to further advance the library are welcomed.
Contributions may not only be in the form of model development, but also
through model use, model testing and validation,
requirements definition or providing feedback regarding the model applicability
to solve specific problems.
</p>
</html>"),
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
        graphics={Bitmap(extent={{-90,-90},{90,90}},
        fileName="modelica://IBPSA/Resources/Images/IBPSA-logo.png")}));
end IBPSA;
