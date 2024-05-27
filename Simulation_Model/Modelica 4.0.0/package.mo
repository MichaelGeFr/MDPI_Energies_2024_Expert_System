package Modelica "Modelica Standard Library - Version 4.0.0"
extends Modelica.Icons.Package;

  package Blocks "Library of basic input/output control blocks (continuous, discrete, logical, table blocks)"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;

    package Continuous "Library of continuous control blocks with internal states"
      import Modelica.Blocks.Interfaces;
      extends Modelica.Icons.Package;

      block Integrator "Output the integral of the input signal with optional reset"
        import Modelica.Blocks.Types.Init;
        parameter Real k(unit="1")=1 "Integrator gain";
        parameter Boolean use_reset = false "= true, if reset port enabled"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true));
        parameter Boolean use_set = false "= true, if set port enabled and used as reinitialization value when reset"
          annotation(Dialog(enable=use_reset), Evaluate=true, HideResult=true, choices(checkBox=true));

        /* InitialState is the default, because it was the default in Modelica 2.2
     and therefore this setting is backward compatible
  */
        parameter Init initType=Init.InitialState
          "Type of initialization (1: no init, 2: steady state, 3,4: initial output)" annotation(Evaluate=true,
            Dialog(group="Initialization"));
        parameter Real y_start=0 "Initial or guess value of output (= state)"
          annotation (Dialog(group="Initialization"));
        extends Interfaces.SISO(y(start=y_start));
        Modelica.Blocks.Interfaces.BooleanInput reset if use_reset "Optional connector of reset signal" annotation(Placement(
          transformation(
            extent={{-20,-20},{20,20}},
            rotation=90,
            origin={60,-120})));
        Modelica.Blocks.Interfaces.RealInput set if use_reset and use_set "Optional connector of set signal" annotation(Placement(
          transformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={60,120})));
    protected
        Modelica.Blocks.Interfaces.BooleanOutput local_reset annotation(HideResult=true);
        Modelica.Blocks.Interfaces.RealOutput local_set annotation(HideResult=true);

      initial equation
        if initType == Init.SteadyState then
           der(y) = 0;
        elseif initType == Init.InitialState or
               initType == Init.InitialOutput then
          y = y_start;
        end if;
      equation
        if use_reset then
          connect(reset, local_reset);
          if use_set then
            connect(set, local_set);
          else
            local_set = y_start;
          end if;
          when local_reset then
            reinit(y, local_set);
          end when;
        else
          local_reset = false;
          local_set = 0;
        end if;
        der(y) = k*u;
        annotation (
          Documentation(info="<html>
<p>
This blocks computes output <strong>y</strong> as
<em>integral</em> of the input <strong>u</strong> multiplied with
the gain <em>k</em>:
</p>
<blockquote><pre>
    k
y = - u
    s
</pre></blockquote>

<p>
It might be difficult to initialize the integrator in steady state.
This is discussed in the description of package
<a href=\"modelica://Modelica.Blocks.Continuous#info\">Continuous</a>.
</p>

<p>
If the <em>reset</em> port is enabled, then the output <strong>y</strong> is reset to <em>set</em>
or to <em>y_start</em> (if the <em>set</em> port is not enabled), whenever the <em>reset</em>
port has a rising edge.
</p>
</html>"),     Icon(coordinateSystem(
                preserveAspectRatio=true,
                extent={{-100.0,-100.0},{100.0,100.0}}),
              graphics={
                Line(
                  points={{-80.0,78.0},{-80.0,-90.0}},
                  color={192,192,192}),
                Polygon(
                  lineColor={192,192,192},
                  fillColor={192,192,192},
                  fillPattern=FillPattern.Solid,
                  points={{-80.0,90.0},{-88.0,68.0},{-72.0,68.0},{-80.0,90.0}}),
                Line(
                  points={{-90.0,-80.0},{82.0,-80.0}},
                  color={192,192,192}),
                Polygon(
                  lineColor={192,192,192},
                  fillColor={192,192,192},
                  fillPattern=FillPattern.Solid,
                  points={{90.0,-80.0},{68.0,-72.0},{68.0,-88.0},{90.0,-80.0}}),
                Text(
                  textColor={192,192,192},
                  extent={{0.0,-70.0},{60.0,-10.0}},
                  textString="I"),
                Text(
                  extent={{-150.0,-150.0},{150.0,-110.0}},
                  textString="k=%k"),
                Line(
                  points=DynamicSelect({{-80.0,-80.0},{80.0,80.0}}, if use_reset then {{-80.0,-80.0},{60.0,60.0},{60.0,-80.0},{80.0,-60.0}} else {{-80.0,-80.0},{80.0,80.0}}),
                  color={0,0,127}),
                Line(
                  visible=use_reset,
                  points={{60,-100},{60,-80}},
                  color={255,0,255},
                  pattern=LinePattern.Dot),
                Text(
                  visible=use_reset,
                  extent={{-28,-62},{94,-86}},
                  textString="reset")}));
      end Integrator;

      block Derivative "Approximated derivative block"
        import Modelica.Blocks.Types.Init;
        parameter Real k(unit="1")=1 "Gains";
        parameter SI.Time T(min=Modelica.Constants.small) = 0.01
          "Time constants (T>0 required; T=0 is ideal derivative block)";
        parameter Init initType=Init.NoInit
          "Type of initialization (1: no init, 2: steady state, 3: initial state, 4: initial output)"
                                                                                          annotation(Evaluate=true,
            Dialog(group="Initialization"));
        parameter Real x_start=0 "Initial or guess value of state"
          annotation (Dialog(group="Initialization"));
        parameter Real y_start=0 "Initial value of output (= state)"
          annotation(Dialog(enable=initType == Init.InitialOutput, group=
                "Initialization"));
        extends Interfaces.SISO;

        output Real x(start=x_start) "State of block";

    protected
        parameter Boolean zeroGain = abs(k) < Modelica.Constants.eps;
      initial equation
        if initType == Init.SteadyState then
          der(x) = 0;
        elseif initType == Init.InitialState then
          x = x_start;
        elseif initType == Init.InitialOutput then
          if zeroGain then
             x = u;
          else
             y = y_start;
          end if;
        end if;
      equation
        der(x) = if zeroGain then 0 else (u - x)/T;
        y = if zeroGain then 0 else (k/T)*(u - x);
        annotation (
          Documentation(info="<html>
<p>
This blocks defines the transfer function between the
input u and the output y
as <em>approximated derivative</em>:
</p>
<blockquote><pre>
        k * s
y = ------------ * u
       T * s + 1
</pre></blockquote>
<p>
If you would like to be able to change easily between different
transfer functions (FirstOrder, SecondOrder, ... ) by changing
parameters, use the general block <strong>TransferFunction</strong> instead
and model a derivative block with parameters<br>
b = {k,0}, a = {T, 1}.
</p>

<p>
If k=0, the block reduces to y=0.
</p>
</html>"),     Icon(
          coordinateSystem(preserveAspectRatio=true,
              extent={{-100.0,-100.0},{100.0,100.0}}),
            graphics={
          Line(points={{-80.0,78.0},{-80.0,-90.0}},
            color={192,192,192}),
        Polygon(lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid,
          points={{-80.0,90.0},{-88.0,68.0},{-72.0,68.0},{-80.0,90.0}}),
        Line(points={{-90.0,-80.0},{82.0,-80.0}},
          color={192,192,192}),
        Polygon(lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid,
          points={{90.0,-80.0},{68.0,-72.0},{68.0,-88.0},{90.0,-80.0}}),
        Line(origin = {-24.667,-27.333},
          points = {{-55.333,87.333},{-19.333,-40.667},{86.667,-52.667}},
          color = {0,0,127},
          smooth = Smooth.Bezier),
        Text(textColor={192,192,192},
          extent={{-30.0,14.0},{86.0,60.0}},
          textString="DT1"),
        Text(extent={{-150.0,-150.0},{150.0,-110.0}},
          textString="k=%k")}));
      end Derivative;

      block FirstOrder "First order transfer function block (= 1 pole)"
        import Modelica.Blocks.Types.Init;
        parameter Real k(unit="1")=1 "Gain";
        parameter SI.Time T(start=1) "Time Constant";
        parameter Init initType=Init.NoInit
          "Type of initialization (1: no init, 2: steady state, 3/4: initial output)" annotation(Evaluate=true,
            Dialog(group="Initialization"));
        parameter Real y_start=0 "Initial or guess value of output (= state)"
          annotation (Dialog(group="Initialization"));

        extends Interfaces.SISO(y(start=y_start));

      initial equation
        if initType == Init.SteadyState then
          der(y) = 0;
        elseif initType == Init.InitialState or initType == Init.InitialOutput then
          y = y_start;
        end if;
      equation
        der(y) = (k*u - y)/T;
        annotation (
          Documentation(info="<html>
<p>
This blocks defines the transfer function between the input u
and the output y as <em>first order</em> system:
</p>
<blockquote><pre>
          k
y = ------------ * u
       T * s + 1
</pre></blockquote>
<p>
If you would like to be able to change easily between different
transfer functions (FirstOrder, SecondOrder, ... ) by changing
parameters, use the general block <strong>TransferFunction</strong> instead
and model a first order SISO system with parameters<br>
b = {k}, a = {T, 1}.
</p>
<blockquote><pre>
Example:
   parameter: k = 0.3, T = 0.4
   results in:
             0.3
      y = ----------- * u
          0.4 s + 1.0
</pre></blockquote>

</html>"),     Icon(
        coordinateSystem(preserveAspectRatio=true,
            extent={{-100.0,-100.0},{100.0,100.0}}),
          graphics={
        Line(points={{-80.0,78.0},{-80.0,-90.0}},
          color={192,192,192}),
        Polygon(lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid,
          points={{-80.0,90.0},{-88.0,68.0},{-72.0,68.0},{-80.0,90.0}}),
        Line(points={{-90.0,-80.0},{82.0,-80.0}},
          color={192,192,192}),
        Polygon(lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid,
          points={{90.0,-80.0},{68.0,-72.0},{68.0,-88.0},{90.0,-80.0}}),
        Line(origin = {-26.667,6.667},
            points = {{106.667,43.333},{-13.333,29.333},{-53.333,-86.667}},
            color = {0,0,127},
            smooth = Smooth.Bezier),
        Text(textColor={192,192,192},
          extent={{0.0,-60.0},{60.0,0.0}},
          textString="PT1"),
        Text(extent={{-150.0,-150.0},{150.0,-110.0}},
          textString="T=%T")}));
      end FirstOrder;

      block LimPID
        "P, PI, PD, and PID controller with limited output, anti-windup compensation, setpoint weighting and optional feed-forward"
        import Modelica.Blocks.Types.Init;
        import Modelica.Blocks.Types.SimpleController;
        extends Modelica.Blocks.Interfaces.SVcontrol;
        output Real controlError = u_s - u_m
          "Control error (set point - measurement)";
        parameter .Modelica.Blocks.Types.SimpleController controllerType=
               .Modelica.Blocks.Types.SimpleController.PID "Type of controller";
        parameter Real k(min=0, unit="1") = 1 "Gain of controller";
        parameter SI.Time Ti(min=Modelica.Constants.small)=0.5
          "Time constant of Integrator block" annotation (Dialog(enable=
                controllerType == .Modelica.Blocks.Types.SimpleController.PI or
                controllerType == .Modelica.Blocks.Types.SimpleController.PID));
        parameter SI.Time Td(min=0)=0.1
          "Time constant of Derivative block" annotation (Dialog(enable=
                controllerType == .Modelica.Blocks.Types.SimpleController.PD or
                controllerType == .Modelica.Blocks.Types.SimpleController.PID));
        parameter Real yMax(start=1) "Upper limit of output";
        parameter Real yMin=-yMax "Lower limit of output";
        parameter Real wp(min=0) = 1
          "Set-point weight for Proportional block (0..1)";
        parameter Real wd(min=0) = 0 "Set-point weight for Derivative block (0..1)"
           annotation(Dialog(enable=controllerType==.Modelica.Blocks.Types.SimpleController.PD or
                                      controllerType==.Modelica.Blocks.Types.SimpleController.PID));
        parameter Real Ni(min=100*Modelica.Constants.eps) = 0.9
          "Ni*Ti is time constant of anti-windup compensation"
           annotation(Dialog(enable=controllerType==.Modelica.Blocks.Types.SimpleController.PI or
                                    controllerType==.Modelica.Blocks.Types.SimpleController.PID));
        parameter Real Nd(min=100*Modelica.Constants.eps) = 10
          "The higher Nd, the more ideal the derivative block"
           annotation(Dialog(enable=controllerType==.Modelica.Blocks.Types.SimpleController.PD or
                                      controllerType==.Modelica.Blocks.Types.SimpleController.PID));
        parameter Boolean withFeedForward=false "Use feed-forward input?"
          annotation(Evaluate=true, choices(checkBox=true));
        parameter Real kFF=1 "Gain of feed-forward input"
          annotation(Dialog(enable=withFeedForward));
        parameter Init initType = Init.InitialState
          "Type of initialization (1: no init, 2: steady state, 3: initial state, 4: initial output)"
          annotation(Evaluate=true, Dialog(group="Initialization"));
        parameter Real xi_start=0
          "Initial or guess value for integrator output (= integrator state)"
          annotation (Dialog(group="Initialization",
                      enable=controllerType==.Modelica.Blocks.Types.SimpleController.PI or
                             controllerType==.Modelica.Blocks.Types.SimpleController.PID));
        parameter Real xd_start=0
          "Initial or guess value for state of derivative block"
          annotation (Dialog(group="Initialization",
                               enable=controllerType==.Modelica.Blocks.Types.SimpleController.PD or
                                      controllerType==.Modelica.Blocks.Types.SimpleController.PID));
        parameter Real y_start=0 "Initial value of output"
          annotation(Dialog(enable=initType == Init.InitialOutput, group=
                "Initialization"));
        parameter Modelica.Blocks.Types.LimiterHomotopy homotopyType = Modelica.Blocks.Types.LimiterHomotopy.Linear
          "Simplified model for homotopy-based initialization"
          annotation (Evaluate=true, Dialog(group="Initialization"));
        parameter Boolean strict=false "= true, if strict limits with noEvent(..)"
          annotation (Evaluate=true, choices(checkBox=true), Dialog(tab="Advanced"));
        constant SI.Time unitTime=1 annotation (HideResult=true);
        Modelica.Blocks.Interfaces.RealInput u_ff if withFeedForward
          "Optional connector of feed-forward input signal"
         annotation (Placement(
              transformation(
              origin={60,-120},
              extent={{20,-20},{-20,20}},
              rotation=270)));
        Modelica.Blocks.Math.Add addP(k1=wp, k2=-1)
          annotation (Placement(transformation(extent={{-80,40},{-60,60}})));
        Modelica.Blocks.Math.Add addD(k1=wd, k2=-1) if with_D
          annotation (Placement(transformation(extent={{-80,-10},{-60,10}})));
        Modelica.Blocks.Math.Gain P(k=1)
          annotation (Placement(transformation(extent={{-50,40},{-30,60}})));
        Modelica.Blocks.Continuous.Integrator I(
          k=unitTime/Ti,
          y_start=xi_start,
          initType=if initType == Init.SteadyState then Init.SteadyState else if
              initType == Init.InitialState
               then Init.InitialState else Init.NoInit) if with_I
          annotation (Placement(transformation(extent={{-50,-60},{-30,-40}})));
        Modelica.Blocks.Continuous.Derivative D(
          k=Td/unitTime,
          T=max([Td/Nd,1.e-14]),
          x_start=xd_start,
          initType=if initType == Init.SteadyState or initType == Init.InitialOutput
               then Init.SteadyState else if initType == Init.InitialState then
              Init.InitialState else Init.NoInit) if with_D
          annotation (Placement(transformation(extent={{-50,-10},{-30,10}})));
        Modelica.Blocks.Math.Gain gainPID(k=k)
          annotation (Placement(transformation(extent={{20,-10},{40,10}})));
        Modelica.Blocks.Math.Add3 addPID
          annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
        Modelica.Blocks.Math.Add3 addI(k2=-1) if with_I
          annotation (Placement(transformation(extent={{-80,-60},{-60,-40}})));
        Modelica.Blocks.Math.Add addSat(k1=+1, k2=-1) if with_I annotation (Placement(
              transformation(
              origin={80,-50},
              extent={{-10,-10},{10,10}},
              rotation=270)));
        Modelica.Blocks.Math.Gain gainTrack(k=1/(k*Ni)) if with_I
          annotation (Placement(transformation(extent={{0,-80},{-20,-60}})));
        Modelica.Blocks.Nonlinear.Limiter limiter(
          uMax=yMax,
          uMin=yMin,
          strict=strict,
          homotopyType=homotopyType)
          annotation (Placement(transformation(extent={{70,-10},{90,10}})));
    protected
        parameter Boolean with_I = controllerType==SimpleController.PI or
                                   controllerType==SimpleController.PID annotation(Evaluate=true, HideResult=true);
        parameter Boolean with_D = controllerType==SimpleController.PD or
                                   controllerType==SimpleController.PID annotation(Evaluate=true, HideResult=true);
    public
        Modelica.Blocks.Sources.Constant Dzero(k=0) if not with_D
          annotation (Placement(transformation(extent={{-40,20},{-30,30}})));
        Modelica.Blocks.Sources.Constant Izero(k=0) if not with_I
          annotation (Placement(transformation(extent={{0,-55},{-10,-45}})));
        Modelica.Blocks.Sources.Constant FFzero(k=0) if not withFeedForward
          annotation (Placement(transformation(extent={{30,-35},{40,-25}})));
        Modelica.Blocks.Math.Add addFF(k1=1, k2=kFF)
          annotation (Placement(transformation(extent={{48,-6},{60,6}})));
      initial equation
        if initType==Init.InitialOutput then
          gainPID.y = y_start;
        end if;
      equation
        if initType == Init.InitialOutput and (y_start < yMin or y_start > yMax) then
            Modelica.Utilities.Streams.error("LimPID: Start value y_start (=" + String(y_start) +
               ") is outside of the limits of yMin (=" + String(yMin) +") and yMax (=" + String(yMax) + ")");
        end if;

        connect(u_s, addP.u1) annotation (Line(points={{-120,0},{-96,0},{-96,56},{
                -82,56}}, color={0,0,127}));
        connect(u_s, addD.u1) annotation (Line(points={{-120,0},{-96,0},{-96,6},{
                -82,6}}, color={0,0,127}));
        connect(u_s, addI.u1) annotation (Line(points={{-120,0},{-96,0},{-96,-42},{
                -82,-42}}, color={0,0,127}));
        connect(addP.y, P.u) annotation (Line(points={{-59,50},{-52,50}}, color={0,
                0,127}));
        connect(addD.y, D.u)
          annotation (Line(points={{-59,0},{-52,0}}, color={0,0,127}));
        connect(addI.y, I.u) annotation (Line(points={{-59,-50},{-52,-50}}, color={
                0,0,127}));
        connect(P.y, addPID.u1) annotation (Line(points={{-29,50},{-20,50},{-20,8},{-12,
                8}},     color={0,0,127}));
        connect(D.y, addPID.u2)
          annotation (Line(points={{-29,0},{-12,0}},color={0,0,127}));
        connect(I.y, addPID.u3) annotation (Line(points={{-29,-50},{-20,-50},{-20,-8},
                {-12,-8}},    color={0,0,127}));
        connect(limiter.y, addSat.u1) annotation (Line(points={{91,0},{94,0},{94,
                -20},{86,-20},{86,-38}}, color={0,0,127}));
        connect(limiter.y, y)
          annotation (Line(points={{91,0},{110,0}}, color={0,0,127}));
        connect(addSat.y, gainTrack.u) annotation (Line(points={{80,-61},{80,-70},{2,-70}},
                          color={0,0,127}));
        connect(gainTrack.y, addI.u3) annotation (Line(points={{-21,-70},{-88,-70},{-88,
                -58},{-82,-58}},     color={0,0,127}));
        connect(u_m, addP.u2) annotation (Line(points={{0,-120},{0,-92},{-92,-92},{-92,44},{-82,44}}, color={0,0,127}));
        connect(u_m, addD.u2) annotation (Line(points={{0,-120},{0,-92},{-92,-92},{-92,-6},{-82,-6}}, color={0,0,127}));
        connect(u_m, addI.u2) annotation (Line(points={{0,-120},{0,-92},{-92,-92},{-92,-50},{-82,-50}}, color={0,0,127}));
        connect(Dzero.y, addPID.u2) annotation (Line(points={{-29.5,25},{-24,25},{-24,
                0},{-12,0}},    color={0,0,127}));
        connect(Izero.y, addPID.u3) annotation (Line(points={{-10.5,-50},{-20,-50},{-20,
                -8},{-12,-8}},    color={0,0,127}));
        connect(addPID.y, gainPID.u)
          annotation (Line(points={{11,0},{18,0}}, color={0,0,127}));
        connect(addFF.y, limiter.u)
          annotation (Line(points={{60.6,0},{68,0}}, color={0,0,127}));
        connect(gainPID.y, addFF.u1) annotation (Line(points={{41,0},{44,0},{44,3.6},
                {46.8,3.6}},color={0,0,127}));
        connect(FFzero.y, addFF.u2) annotation (Line(points={{40.5,-30},{44,-30},{44,
                -3.6},{46.8,-3.6}},
                              color={0,0,127}));
        connect(addFF.u2, u_ff) annotation (Line(points={{46.8,-3.6},{44,-3.6},{44,
                -92},{60,-92},{60,-120}},
                                     color={0,0,127}));
        connect(addFF.y, addSat.u2) annotation (Line(points={{60.6,0},{64,0},{64,-20},
                {74,-20},{74,-38}}, color={0,0,127}));
        annotation (defaultComponentName="PID",
          Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{-80,78},{-80,-90}}, color={192,192,192}),
              Polygon(
                points={{-80,90},{-88,68},{-72,68},{-80,90}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-90,-80},{82,-80}}, color={192,192,192}),
              Polygon(
                points={{90,-80},{68,-72},{68,-88},{90,-80}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-80,-80},{-80,-20},{30,60},{80,60}}, color={0,0,127}),
              Text(
                extent={{-20,-20},{80,-60}},
                textColor={192,192,192},
                textString="%controllerType"),
              Line(
                visible=strict,
                points={{30,60},{81,60}},
                color={255,0,0})}),
          Diagram(graphics={Text(
                  extent={{79,-112},{129,-102}},
                  textColor={0,0,255},
                textString=" (feed-forward)")}),
          Documentation(info="<html>
<p>
Via parameter <strong>controllerType</strong> either <strong>P</strong>, <strong>PI</strong>, <strong>PD</strong>,
or <strong>PID</strong> can be selected. If, e.g., PI is selected, all components belonging to the
D-part are removed from the block (via conditional declarations).
The example model
<a href=\"modelica://Modelica.Blocks.Examples.PID_Controller\">Modelica.Blocks.Examples.PID_Controller</a>
demonstrates the usage of this controller.
Several practical aspects of PID controller design are incorporated
according to chapter 3 of the book:
</p>

<dl>
<dt>&Aring;str&ouml;m K.J., and H&auml;gglund T.:</dt>
<dd> <strong>PID Controllers: Theory, Design, and Tuning</strong>.
     Instrument Society of America, 2nd edition, 1995.
</dd>
</dl>

<p>
Besides the additive <strong>proportional, integral</strong> and <strong>derivative</strong>
part of this controller, the following features are present:
</p>
<ul>
<li> The output of this controller is limited. If the controller is
     in its limits, anti-windup compensation is activated to drive
     the integrator state to zero.</li>
<li> The high-frequency gain of the derivative part is limited
     to avoid excessive amplification of measurement noise.</li>
<li> Setpoint weighting is present, which allows to weight
     the setpoint in the proportional and the derivative part
     independently from the measurement. The controller will respond
     to load disturbances and measurement noise independently of this setting
     (parameters wp, wd). However, setpoint changes will depend on this
     setting. For example, it is useful to set the setpoint weight wd
     for the derivative part to zero, if steps may occur in the
     setpoint signal.</li>
<li> Optional feed-forward. It is possible to add a feed-forward signal.
     The feed-forward signal is added before limitation.</li>
</ul>

<p>
The parameters of the controller can be manually adjusted by performing
simulations of the closed loop system (= controller + plant connected
together) and using the following strategy:
</p>

<ol>
<li> Set very large limits, e.g., yMax = Modelica.Constants.inf</li>
<li> Select a <strong>P</strong>-controller and manually enlarge parameter <strong>k</strong>
     (the total gain of the controller) until the closed-loop response
     cannot be improved any more.</li>
<li> Select a <strong>PI</strong>-controller and manually adjust parameters
     <strong>k</strong> and <strong>Ti</strong> (the time constant of the integrator).
     The first value of Ti can be selected, such that it is in the
     order of the time constant of the oscillations occurring with
     the P-controller. If, e.g., vibrations in the order of T=10 ms
     occur in the previous step, start with Ti=0.01 s.</li>
<li> If you want to make the reaction of the control loop faster
     (but probably less robust against disturbances and measurement noise)
     select a <strong>PID</strong>-Controller and manually adjust parameters
     <strong>k</strong>, <strong>Ti</strong>, <strong>Td</strong> (time constant of derivative block).</li>
<li> Set the limits yMax and yMin according to your specification.</li>
<li> Perform simulations such that the output of the PID controller
     goes in its limits. Tune <strong>Ni</strong> (Ni*Ti is the time constant of
     the anti-windup compensation) such that the input to the limiter
     block (= limiter.u) goes quickly enough back to its limits.
     If Ni is decreased, this happens faster. If Ni=infinity, the
     anti-windup compensation is switched off and the controller works bad.</li>
</ol>

<p>
<strong>Initialization</strong>
</p>

<p>
This block can be initialized in different
ways controlled by parameter <strong>initType</strong>. The possible
values of initType are defined in
<a href=\"modelica://Modelica.Blocks.Types.Init\">Modelica.Blocks.Types.Init</a>.
</p>

<p>
Based on the setting of initType, the integrator (I) and derivative (D)
blocks inside the PID controller are initialized according to the following table:
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>initType</strong></td>
      <td><strong>I.initType</strong></td>
      <td><strong>D.initType</strong></td></tr>

  <tr><td><strong>NoInit</strong></td>
      <td>NoInit</td>
      <td>NoInit</td></tr>

  <tr><td><strong>SteadyState</strong></td>
      <td>SteadyState</td>
      <td>SteadyState</td></tr>

  <tr><td><strong>InitialState</strong></td>
      <td>InitialState</td>
      <td>InitialState</td></tr>

  <tr><td><strong>InitialOutput</strong><br>
          and initial equation: y = y_start</td>
      <td>NoInit</td>
      <td>SteadyState</td></tr>
</table>

<p>
In many cases, the most useful initial condition is
<strong>SteadyState</strong> because initial transients are then no longer
present. If initType = Init.SteadyState, then in some
cases difficulties might occur. The reason is the
equation of the integrator:
</p>

<blockquote><pre>
<strong>der</strong>(y) = k*u;
</pre></blockquote>

<p>
The steady state equation \"der(x)=0\" leads to the condition that the input u to the
integrator is zero. If the input u is already (directly or indirectly) defined
by another initial condition, then the initialization problem is <strong>singular</strong>
(has none or infinitely many solutions). This situation occurs often
for mechanical systems, where, e.g., u = desiredSpeed - measuredSpeed and
since speed is both a state and a derivative, it is natural to
initialize it with zero. As sketched this is, however, not possible.
The solution is to not initialize u_m or the variable that is used
to compute u_m by an algebraic equation.
</p>

<p>
When initializing in steady-state, homotopy-based initialization can help the convergence of the solver,
by using a simplified model a the beginning of the solution process. Different options are available.
</p>

<ul>
<li><strong>homotopyType=Linear</strong> (default): the limitations are removed from the simplified model,
making it linear. Use this if you know that the controller will not be saturated at steady state.</li>
<li><strong>homotopyType=UpperLimit</strong>: if it is known a priori the controller will be stuck at the upper
limit yMax, this option assumes y = yMax as a simplified model.</li>
<li><strong>homotopyType=LowerLimit</strong>: if it is known a priori the controller will be stuck at the lower
limit yMin, this option assumes y = yMin as a simplified model.</li>
<li><strong>homotopyType=NoHomotopy</strong>: this option does not apply any simplification and keeps the
limiter active throughout the homotopy transformation. Use this if it is unknown whether the controller
is saturated or not at initialization and if the limitations on the output must be enforced throughout
the entire homotopy transformation.</li>
</ul>
</html>"));
      end LimPID;
      annotation (
        Documentation(info="<html>
<p>
This package contains basic <strong>continuous</strong> input/output blocks
described by differential equations.
</p>

<p>
All blocks of this package can be initialized in different
ways controlled by parameter <strong>initType</strong>. The possible
values of initType are defined in
<a href=\"modelica://Modelica.Blocks.Types.Init\">Modelica.Blocks.Types.Init</a>:
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>Name</strong></td>
      <td><strong>Description</strong></td></tr>

  <tr><td><strong>Init.NoInit</strong></td>
      <td>no initialization (start values are used as guess values with fixed=false)</td></tr>

  <tr><td><strong>Init.SteadyState</strong></td>
      <td>steady state initialization (derivatives of states are zero)</td></tr>

  <tr><td><strong>Init.InitialState</strong></td>
      <td>Initialization with initial states</td></tr>

  <tr><td><strong>Init.InitialOutput</strong></td>
      <td>Initialization with initial outputs (and steady state of the states if possible)</td></tr>
</table>

<p>
For backward compatibility reasons the default of all blocks is
<strong>Init.NoInit</strong>, with the exception of Integrator and LimIntegrator
where the default is <strong>Init.InitialState</strong> (this was the initialization
defined in version 2.2 of the Modelica standard library).
</p>

<p>
In many cases, the most useful initial condition is
<strong>Init.SteadyState</strong> because initial transients are then no longer
present. The drawback is that in combination with a non-linear
plant, non-linear algebraic equations occur that might be
difficult to solve if appropriate guess values for the
iteration variables are not provided (i.e., start values with fixed=false).
However, it is often already useful to just initialize
the linear blocks from the Continuous blocks library in SteadyState.
This is uncritical, because only linear algebraic equations occur.
If Init.NoInit is set, then the start values for the states are
interpreted as <strong>guess</strong> values and are propagated to the
states with fixed=<strong>false</strong>.
</p>

<p>
Note, initialization with Init.SteadyState is usually difficult
for a block that contains an integrator
(Integrator, LimIntegrator, PI, PID, LimPID).
This is due to the basic equation of an integrator:
</p>

<blockquote><pre>
<strong>initial equation</strong>
   <strong>der</strong>(y) = 0;   // Init.SteadyState
<strong>equation</strong>
   <strong>der</strong>(y) = k*u;
</pre></blockquote>

<p>
The steady state equation leads to the condition that the input to the
integrator is zero. If the input u is already (directly or indirectly) defined
by another initial condition, then the initialization problem is <strong>singular</strong>
(has none or infinitely many solutions). This situation occurs often
for mechanical systems, where, e.g., u = desiredSpeed - measuredSpeed and
since speed is both a state and a derivative, it is always defined by
Init.InitialState or Init.SteadyState initialization.
</p>

<p>
In such a case, <strong>Init.NoInit</strong> has to be selected for the integrator
and an additional initial equation has to be added to the system
to which the integrator is connected. E.g., useful initial conditions
for a 1-dim. rotational inertia controlled by a PI controller are that
<strong>angle</strong>, <strong>speed</strong>, and <strong>acceleration</strong> of the inertia are zero.
</p>

</html>"),     Icon(graphics={Line(
              origin={0.061,4.184},
              points={{81.939,36.056},{65.362,36.056},{14.39,-26.199},{-29.966,
                  113.485},{-65.374,-61.217},{-78.061,-78.184}},
              color={95,95,95},
              smooth=Smooth.Bezier)}));
    end Continuous;

    package Interfaces "Library of connectors and partial models for input/output blocks"
      extends Modelica.Icons.InterfacesPackage;

      connector RealInput = input Real "'input Real' as connector" annotation (
        defaultComponentName="u",
        Icon(graphics={
          Polygon(
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid,
            points={{-100.0,100.0},{100.0,0.0},{-100.0,-100.0}})},
          coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}},
            preserveAspectRatio=true,
            initialScale=0.2)),
        Diagram(
          coordinateSystem(preserveAspectRatio=true,
            initialScale=0.2,
            extent={{-100.0,-100.0},{100.0,100.0}}),
            graphics={
          Polygon(
            lineColor={0,0,127},
            fillColor={0,0,127},
            fillPattern=FillPattern.Solid,
            points={{0.0,50.0},{100.0,0.0},{0.0,-50.0},{0.0,50.0}}),
          Text(
            textColor={0,0,127},
            extent={{-10.0,60.0},{-10.0,85.0}},
            textString="%name")}),
        Documentation(info="<html>
<p>
Connector with one input signal of type Real.
</p>
</html>"));

      connector RealOutput = output Real "'output Real' as connector" annotation (
        defaultComponentName="y",
        Icon(
          coordinateSystem(preserveAspectRatio=true,
            extent={{-100.0,-100.0},{100.0,100.0}}),
            graphics={
          Polygon(
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            points={{-100.0,100.0},{100.0,0.0},{-100.0,-100.0}})}),
        Diagram(
          coordinateSystem(preserveAspectRatio=true,
            extent={{-100.0,-100.0},{100.0,100.0}}),
            graphics={
          Polygon(
            lineColor={0,0,127},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid,
            points={{-100.0,50.0},{0.0,0.0},{-100.0,-50.0}}),
          Text(
            textColor={0,0,127},
            extent={{30.0,60.0},{30.0,110.0}},
            textString="%name")}),
        Documentation(info="<html>
<p>
Connector with one output signal of type Real.
</p>
</html>"));

      connector BooleanInput = input Boolean "'input Boolean' as connector"
        annotation (
        defaultComponentName="u",
        Icon(graphics={Polygon(
              points={{-100,100},{100,0},{-100,-100},{-100,100}},
              lineColor={255,0,255},
              fillColor={255,0,255},
              fillPattern=FillPattern.Solid)}, coordinateSystem(
            extent={{-100,-100},{100,100}},
            preserveAspectRatio=true,
            initialScale=0.2)),
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            initialScale=0.2,
            extent={{-100,-100},{100,100}}), graphics={Polygon(
              points={{0,50},{100,0},{0,-50},{0,50}},
              lineColor={255,0,255},
              fillColor={255,0,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{-10,85},{-10,60}},
              textColor={255,0,255},
              textString="%name")}),
        Documentation(info="<html>
<p>
Connector with one input signal of type Boolean.
</p>
</html>"));

      connector BooleanOutput = output Boolean "'output Boolean' as connector"
        annotation (
        defaultComponentName="y",
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={Polygon(
              points={{-100,100},{100,0},{-100,-100},{-100,100}},
              lineColor={255,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}}), graphics={Polygon(
              points={{-100,50},{0,0},{-100,-50},{-100,50}},
              lineColor={255,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid), Text(
              extent={{30,110},{30,60}},
              textColor={255,0,255},
              textString="%name")}),
        Documentation(info="<html>
<p>
Connector with one output signal of type Boolean.
</p>
</html>"));

      connector RealVectorInput = input Real
        "Real input connector used for vector of connectors" annotation (
        defaultComponentName="u",
        Icon(graphics={Ellipse(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid)}, coordinateSystem(
            extent={{-100,-100},{100,100}},
            preserveAspectRatio=true,
            initialScale=0.2)),
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            initialScale=0.2,
            extent={{-100,-100},{100,100}}), graphics={Text(
              extent={{-10,85},{-10,60}},
              textColor={0,0,127},
              textString="%name"), Ellipse(
              extent={{-50,50},{50,-50}},
              lineColor={0,0,127},
              fillColor={0,0,127},
              fillPattern=FillPattern.Solid)}),
        Documentation(info="<html>
<p>
Real input connector that is used for a vector of connectors,
for example <a href=\"modelica://Modelica.Blocks.Interfaces.PartialRealMISO\">PartialRealMISO</a>,
and has therefore a different icon as RealInput connector.
</p>
</html>"));

      partial block SO "Single Output continuous control block"
        extends Modelica.Blocks.Icons.Block;

        RealOutput y "Connector of Real output signal" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));
        annotation (Documentation(info="<html>
<p>
Block has one continuous Real output signal.
</p>
</html>"));

      end SO;

      partial block MO "Multiple Output continuous control block"
        extends Modelica.Blocks.Icons.Block;

        parameter Integer nout(min=1) = 1 "Number of outputs";
        RealOutput y[nout] "Connector of Real output signals" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));
        annotation (Documentation(info="<html>
<p>
Block has one continuous Real output signal vector.
</p>
</html>"));

      end MO;

      partial block SISO "Single Input Single Output continuous control block"
        extends Modelica.Blocks.Icons.Block;

        RealInput u "Connector of Real input signal" annotation (Placement(
              transformation(extent={{-140,-20},{-100,20}})));
        RealOutput y "Connector of Real output signal" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));
        annotation (Documentation(info="<html>
<p>
Block has one continuous Real input and one continuous Real output signal.
</p>
</html>"));
      end SISO;

      partial block SI2SO
        "2 Single Input / 1 Single Output continuous control block"
        extends Modelica.Blocks.Icons.Block;

        RealInput u1 "Connector of Real input signal 1" annotation (Placement(
              transformation(extent={{-140,40},{-100,80}})));
        RealInput u2 "Connector of Real input signal 2" annotation (Placement(
              transformation(extent={{-140,-80},{-100,-40}})));
        RealOutput y "Connector of Real output signal" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));

        annotation (Documentation(info="<html>
<p>
Block has two continuous Real input signals u1 and u2 and one
continuous Real output signal y.
</p>
</html>"));

      end SI2SO;

      partial block PartialRealMISO
        "Partial block with a RealVectorInput and a RealOutput signal"

        parameter Integer significantDigits(min=1) = 3
          "Number of significant digits to be shown in dynamic diagram layer for y"
          annotation (Dialog(tab="Advanced"));
        parameter Integer nu(min=0) = 0 "Number of input connections"
          annotation (Dialog(connectorSizing=true), HideResult=true);
        Modelica.Blocks.Interfaces.RealVectorInput u[nu]
          annotation (Placement(transformation(extent={{-120,70},{-80,-70}})));
        Modelica.Blocks.Interfaces.RealOutput y
          annotation (Placement(transformation(extent={{100,-17},{134,17}})));
        annotation (Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              initialScale=0.06), graphics={
              Text(
                extent={{110,-50},{300,-70}},
                textString=DynamicSelect(" ", String(y, significantDigits=
                    significantDigits))),
              Text(
                extent={{-250,170},{250,110}},
                textString="%name",
                textColor={0,0,255}),
              Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={255,137,0},
                fillColor={255,255,255},
                borderPattern=BorderPattern.Raised,
                fillPattern=FillPattern.Solid)}));
      end PartialRealMISO;

      partial block SVcontrol "Single-Variable continuous controller"
        extends Modelica.Blocks.Icons.Block;

        RealInput u_s "Connector of setpoint input signal" annotation (Placement(
              transformation(extent={{-140,-20},{-100,20}})));
        RealInput u_m "Connector of measurement input signal" annotation (Placement(
              transformation(
              origin={0,-120},
              extent={{20,-20},{-20,20}},
              rotation=270)));
        RealOutput y "Connector of actuator output signal" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));
        annotation (Documentation(info="<html>
<p>
Block has two continuous Real input signals and one
continuous Real output signal. The block is designed
to be used as base class for a corresponding controller.
</p>
</html>"));
      end SVcontrol;

      partial block partialBooleanSI "Partial block with 1 input Boolean signal"
        extends Modelica.Blocks.Icons.PartialBooleanBlock;

        Blocks.Interfaces.BooleanInput u "Connector of Boolean input signal"
          annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics={Ellipse(
                extent={{-71,7},{-85,-7}},
                lineColor=DynamicSelect({235,235,235}, if u then {0,255,0} else {235,235,235}),
                fillColor=DynamicSelect({235,235,235}, if u then {0,255,0} else {235,235,235}),
                fillPattern=FillPattern.Solid)}), Documentation(info="<html>
<p>
Block has one continuous Boolean input signal
with a 3D icon (e.g., used in Blocks.Logical library).
</p>
</html>"));

      end partialBooleanSI;

      partial block partialBooleanSO "Partial block with 1 output Boolean signal"

        Blocks.Interfaces.BooleanOutput y "Connector of Boolean output signal"
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));
        extends Modelica.Blocks.Icons.PartialBooleanBlock;

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics={Ellipse(
                extent={{71,7},{85,-7}},
                lineColor=DynamicSelect({235,235,235}, if y then {0,255,0} else {235,235,235}),
                fillColor=DynamicSelect({235,235,235}, if y then {0,255,0} else {235,235,235}),
                fillPattern=FillPattern.Solid)}), Documentation(info="<html>
<p>
Block has one continuous Boolean output signal
with a 3D icon (e.g., used in Blocks.Logical library).
</p>
</html>"));

      end partialBooleanSO;

      partial block PartialBooleanSISO_small
        "Partial block with a BooleanInput and a BooleanOutput signal and a small block icon"

        Modelica.Blocks.Interfaces.BooleanInput u "Boolean input signal"
          annotation (Placement(transformation(extent={{-180,-40},{-100,40}})));
        Modelica.Blocks.Interfaces.BooleanOutput y "Boolean output signal"
          annotation (Placement(transformation(extent={{100,-20},{140,20}})));
        annotation (Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              initialScale=0.04), graphics={
              Text(
                extent={{-300,200},{300,120}},
                textString="%name",
                textColor={0,0,255}),
              Rectangle(
                extent={{-100,100},{100,-100}},
                fillColor={215,215,215},
                fillPattern=FillPattern.Solid,
                borderPattern=BorderPattern.Raised),
              Ellipse(
                extent={{60,10},{80,-10}},
                lineColor=DynamicSelect({235,235,235}, if y then {0,255,0} else {235,235,235}),
                fillColor=DynamicSelect({235,235,235}, if y then {0,255,0} else {235,235,235}),
                fillPattern=FillPattern.Solid)}));
      end PartialBooleanSISO_small;
      annotation (Documentation(info="<html>
<p>
This package contains interface definitions for
<strong>continuous</strong> input/output blocks with Real,
Integer and Boolean signals. Furthermore, it contains
partial models for continuous and discrete blocks.
</p>

</html>",     revisions="<html>
<ul>
<li><em>June 28, 2019</em>
       by Thomas Beutlich:<br>
       Removed obsolete blocks.</li>
<li><em>Oct. 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       Added several new interfaces.</li>
<li><em>Oct. 24, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       RealInputSignal renamed to RealInput. RealOutputSignal renamed to
       output RealOutput. GraphBlock renamed to BlockIcon. SISOreal renamed to
       SISO. SOreal renamed to SO. I2SOreal renamed to M2SO.
       SignalGenerator renamed to SignalSource. Introduced the following
       new models: MIMO, MIMOs, SVcontrol, MVcontrol, DiscreteBlockIcon,
       DiscreteBlock, DiscreteSISO, DiscreteMIMO, DiscreteMIMOs,
       BooleanBlockIcon, BooleanSISO, BooleanSignalSource, MI2BooleanMOs.</li>
<li><em>June 30, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized a first version, based on an existing Dymola library
       of Dieter Moormann and Hilding Elmqvist.</li>
</ul>
</html>"));
    end Interfaces;

    package Math "Library of Real mathematical functions as input/output blocks"
      import Modelica.Blocks.Interfaces;
      extends Modelica.Icons.Package;

      block Gain "Output the product of a gain value with the input signal"

        parameter Real k(start=1, unit="1")
          "Gain value multiplied with input signal";
    public
        Interfaces.RealInput u "Input signal connector" annotation (Placement(
              transformation(extent={{-140,-20},{-100,20}})));
        Interfaces.RealOutput y "Output signal connector" annotation (Placement(
              transformation(extent={{100,-10},{120,10}})));

      equation
        y = k*u;
        annotation (
          Documentation(info="<html>
<p>
This block computes output <em>y</em> as
<em>product</em> of gain <em>k</em> with the
input <em>u</em>:
</p>
<blockquote><pre>
y = k * u;
</pre></blockquote>

</html>"),Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                  100}}), graphics={
              Polygon(
                points={{-100,-100},{-100,100},{100,0},{-100,-100}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-150,-140},{150,-100}},
                textString="k=%k"),
              Text(
                extent={{-150,140},{150,100}},
                textString="%name",
                textColor={0,0,255})}));
      end Gain;

      block MultiSum "Sum of Reals: y = k[1]*u[1] + k[2]*u[2] + ... + k[n]*u[n]"
        extends Modelica.Blocks.Interfaces.PartialRealMISO;
        parameter Real k[nu]=fill(1, nu) "Input gains";
      equation
        if size(u, 1) > 0 then
          y = k*u;
        else
          y = 0;
        end if;

        annotation (Icon(graphics={Text(
                extent={{-200,-110},{200,-140}},
                textString="%k"), Text(
                extent={{-72,68},{92,-68}},
                textString="+")}), Documentation(info="<html>
<p>
This blocks computes the scalar Real output \"y\" as sum of the elements of the
Real input signal vector u:
</p>
<blockquote><pre>
y = k[1]*u[1] + k[2]*u[2] + ... k[N]*u[N];
</pre></blockquote>

<p>
The input connector is a vector of Real input signals.
When a connection line is drawn, the dimension of the input
vector is enlarged by one and the connection is automatically
connected to this new free index (thanks to the
connectorSizing annotation).
</p>

<p>
The usage is demonstrated, e.g., in example
<a href=\"modelica://Modelica.Blocks.Examples.RealNetwork1\">Modelica.Blocks.Examples.RealNetwork1</a>.
</p>

<p>
If no connection to the input connector \"u\" is present,
the output is set to zero: y=0.
</p>

</html>"));
      end MultiSum;

      block Add "Output the sum of the two inputs"
        extends Interfaces.SI2SO;

        parameter Real k1=+1 "Gain of input signal 1";
        parameter Real k2=+1 "Gain of input signal 2";

      equation
        y = k1*u1 + k2*u2;
        annotation (
          Documentation(info="<html>
<p>
This blocks computes output <strong>y</strong> as <em>sum</em> of the
two input signals <strong>u1</strong> and <strong>u2</strong>:
</p>
<blockquote><pre>
<strong>y</strong> = k1*<strong>u1</strong> + k2*<strong>u2</strong>;
</pre></blockquote>
<p>
Example:
</p>
<blockquote><pre>
   parameter:   k1= +2, k2= -3

results in the following equations:

   y = 2 * u1 - 3 * u2
</pre></blockquote>

</html>"),Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{-100,60},{-74,24},{-44,24}}, color={0,0,127}),
              Line(points={{-100,-60},{-74,-24},{-44,-24}}, color={0,0,127}),
              Ellipse(lineColor={0,0,127}, extent={{-50,-50},{50,50}}),
              Line(points={{50,0},{100,0}}, color={0,0,127}),
              Text(extent={{-40,40},{40,-40}}, textString="+"),
              Text(extent={{-100,52},{5,92}}, textString="%k1"),
              Text(extent={{-100,-92},{5,-52}}, textString="%k2")}),
          Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                  100,100}}), graphics={         Line(points={{50,0},{100,0}},
                color={0,0,255}),                                        Line(
                points={{50,0},{100,0}}, color={0,0,127})}));
      end Add;

      block Add3 "Output the sum of the three inputs"
        extends Modelica.Blocks.Icons.Block;

        parameter Real k1=+1 "Gain of input signal 1";
        parameter Real k2=+1 "Gain of input signal 2";
        parameter Real k3=+1 "Gain of input signal 3";
        Interfaces.RealInput u1 "Connector of Real input signal 1" annotation (
            Placement(transformation(extent={{-140,60},{-100,100}})));
        Interfaces.RealInput u2 "Connector of Real input signal 2" annotation (
            Placement(transformation(extent={{-140,-20},{-100,20}})));
        Interfaces.RealInput u3 "Connector of Real input signal 3" annotation (
            Placement(transformation(extent={{-140,-100},{-100,-60}})));
        Interfaces.RealOutput y "Connector of Real output signal" annotation (
            Placement(transformation(extent={{100,-10},{120,10}})));

      equation
        y = k1*u1 + k2*u2 + k3*u3;
        annotation (
          Documentation(info="<html>
<p>
This blocks computes output <strong>y</strong> as <em>sum</em> of the
three input signals <strong>u1</strong>, <strong>u2</strong> and <strong>u3</strong>:
</p>
<blockquote><pre>
<strong>y</strong> = k1*<strong>u1</strong> + k2*<strong>u2</strong> + k3*<strong>u3</strong>;
</pre></blockquote>
<p>
Example:
</p>
<blockquote><pre>
   parameter:   k1= +2, k2= -3, k3=1;

results in the following equations:

   y = 2 * u1 - 3 * u2 + u3;
</pre></blockquote>

</html>"),Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                  100}}), graphics={
              Text(
                extent={{-100,50},{5,90}},
                textString="%k1"),
              Text(
                extent={{-100,-20},{5,20}},
                textString="%k2"),
              Text(
                extent={{-100,-50},{5,-90}},
                textString="%k3"),
              Text(
                extent={{10,40},{90,-40}},
                textString="+")}));
      end Add3;

      block Product "Output product of the two inputs"
        extends Interfaces.SI2SO;

      equation
        y = u1*u2;
        annotation (
          Documentation(info="<html>
<p>
This blocks computes the output <strong>y</strong>
as <em>product</em> of the two inputs <strong>u1</strong> and <strong>u2</strong>:
</p>
<blockquote><pre>
y = u1 * u2;
</pre></blockquote>

</html>"),Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{-100,60},{-40,60},{-30,40}}, color={0,0,127}),
              Line(points={{-100,-60},{-40,-60},{-30,-40}}, color={0,0,127}),
              Line(points={{50,0},{100,0}}, color={0,0,127}),
              Line(points={{-30,0},{30,0}}),
              Line(points={{-15,25.99},{15,-25.99}}),
              Line(points={{-15,-25.99},{15,25.99}}),
              Ellipse(lineColor={0,0,127}, extent={{-50,-50},{50,50}})}));
      end Product;

      block BooleanToReal "Convert Boolean to Real signal"
        extends Interfaces.partialBooleanSI;
        parameter Real realTrue=1.0 "Output signal for true Boolean input";
        parameter Real realFalse=0.0 "Output signal for false Boolean input";

        Blocks.Interfaces.RealOutput y "Connector of Real output signal"
          annotation (Placement(transformation(extent={{100,-10},{120,10}})));

      equation
        y = if u then realTrue else realFalse;
        annotation (Documentation(info="<html>
<p>
This block computes the output <strong>y</strong>
as <em>Real equivalent</em> of the Boolean input <strong>u</strong>:
</p>
<blockquote><pre>
y = <strong>if</strong> u <strong>then</strong> realTrue <strong>else</strong> realFalse;
</pre></blockquote>
<p>where <strong>u</strong> is of Boolean and <strong>y</strong> of Real type,
and <strong>realTrue</strong> and <strong>realFalse</strong> are parameters.
</p>
</html>"),     Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics={
              Text(
                extent={{-86,92},{-6,10}},
                textColor={255,0,255},
                textString="B"),
              Polygon(
                points={{-12,-46},{-32,-26},{-32,-36},{-64,-36},{-64,-56},{-32,-56},
                    {-32,-66},{-12,-46}},
                fillColor={0,0,127},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,127}),
              Text(
                extent={{8,-4},{92,-94}},
                textString="R",
                textColor={0,0,127})}));
      end BooleanToReal;

      block RealToBoolean "Convert Real to Boolean signal"

        Blocks.Interfaces.RealInput u "Connector of Real input signal" annotation (
            Placement(transformation(extent={{-140,-20},{-100,20}})));
        extends Interfaces.partialBooleanSO;
        parameter Real threshold=0.5
          "Output signal y is true, if input u >= threshold";

      equation
        y = u >= threshold;
        annotation (Documentation(info="<html>
<p>
This block computes the Boolean output <strong>y</strong>
from the Real input <strong>u</strong> by the equation:
</p>

<blockquote><pre>
y = u &ge; threshold;
</pre></blockquote>

<p>
where <strong>threshold</strong> is a parameter.
</p>
</html>"),     Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                  {100,100}}), graphics={
              Text(
                extent={{-86,92},{-6,10}},
                textColor={0,0,127},
                textString="R"),
              Polygon(
                points={{-12,-46},{-32,-26},{-32,-36},{-64,-36},{-64,-56},{-32,-56},
                    {-32,-66},{-12,-46}},
                lineColor={255,0,255},
                fillColor={255,0,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{8,-4},{92,-94}},
                textColor={255,0,255},
                textString="B")}));
      end RealToBoolean;
      annotation (Documentation(info="<html>
<p>
This package contains basic <strong>mathematical operations</strong>,
such as summation and multiplication, and basic <strong>mathematical
functions</strong>, such as <strong>sqrt</strong> and <strong>sin</strong>, as
input/output blocks. All blocks of this library can be either
connected with continuous blocks or with sampled-data blocks.
</p>
</html>",     revisions="<html>
<ul>
<li><em>August 24, 2016</em>
       by Christian Kral: added WrapAngle</li>
<li><em>October 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       New blocks added: RealToInteger, IntegerToReal, Max, Min, Edge, BooleanChange, IntegerChange.</li>
<li><em>August 7, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized (partly based on an existing Dymola library
       of Dieter Moormann and Hilding Elmqvist).
</li>
</ul>
</html>"),     Icon(graphics={Line(
              points={{-80,-2},{-68.7,32.2},{-61.5,51.1},{-55.1,64.4},{-49.4,72.6},
                  {-43.8,77.1},{-38.2,77.8},{-32.6,74.6},{-26.9,67.7},{-21.3,57.4},
                  {-14.9,42.1},{-6.83,19.2},{10.1,-32.8},{17.3,-52.2},{23.7,-66.2},
                  {29.3,-75.1},{35,-80.4},{40.6,-82},{46.2,-79.6},{51.9,-73.5},{
                  57.5,-63.9},{63.9,-49.2},{72,-26.8},{80,-2}},
              color={95,95,95},
              smooth=Smooth.Bezier)}));
    end Math;

    package MathBoolean "Library of Boolean mathematical functions as input/output blocks"
      extends Modelica.Icons.Package;

      block Not "Logical 'not': y = not u"
        extends Modelica.Blocks.Interfaces.PartialBooleanSISO_small;

      equation
        y = not u;
        annotation (defaultComponentName="not1", Icon(graphics={Text(
                  extent={{-98,40},{42,-40}},
                  textString="not")}),
          Documentation(info="<html>
<p>
The output is <strong>false</strong> if at least one input is <strong>true</strong>, otherwise
the output is <strong>true</strong>.
</p>

<p>
The input connector is a vector of Boolean input signals.
When a connection line is drawn, the dimension of the input
vector is enlarged by one and the connection is automatically
connected to this new free index (thanks to the
connectorSizing annotation).
</p>

<p>
The usage is demonstrated, e.g., in example
<a href=\"modelica://Modelica.Blocks.Examples.BooleanNetwork1\">Modelica.Blocks.Examples.BooleanNetwork1</a>.
</p>
</html>"));
      end Not;
      annotation (Documentation(info="<html>
<p>
This package contains basic <strong>mathematical operations</strong>
on <strong>Boolean</strong> signals.
</p>

<p>
The new features are:
</p>

<ul>
<li> If useful, blocks may have an arbitrary number of inputs (e.g., \"And\" block with 2,3,4,...
     Boolean inputs). This is based on the \"connectorSizing\" annotation which
     allows a tool to conveniently handle vectors of connectors.</li>

<li> The blocks are smaller in size, so that the diagram area is better
     utilized for trivial blocks such as \"And\" or \"Or\".</li>

</ul>

</html>"),     Icon(graphics={Line(points={{-80,-16},{-4,-16},{-4,28},{38,28},{38,
                  -16},{66,-16}}, color={255,0,255})}));
    end MathBoolean;

    package Nonlinear "Library of discontinuous or non-differentiable algebraic control blocks"
      import Modelica.Blocks.Interfaces;
      extends Modelica.Icons.Package;

          block Limiter "Limit the range of a signal"
            parameter Real uMax(start=1) "Upper limits of input signals";
            parameter Real uMin= -uMax "Lower limits of input signals";
            parameter Boolean strict=false "= true, if strict limits with noEvent(..)"
              annotation (Evaluate=true, choices(checkBox=true), Dialog(tab="Advanced"));
            parameter Types.LimiterHomotopy homotopyType = Modelica.Blocks.Types.LimiterHomotopy.Linear "Simplified model for homotopy-based initialization"
              annotation (Evaluate=true, Dialog(group="Initialization"));
            extends Interfaces.SISO;
    protected
            Real simplifiedExpr "Simplified expression for homotopy-based initialization";

          equation
            assert(uMax >= uMin, "Limiter: Limits must be consistent. However, uMax (=" + String(uMax) +
                                 ") < uMin (=" + String(uMin) + ")");
            simplifiedExpr = (if homotopyType == Types.LimiterHomotopy.Linear then u
                              else if homotopyType == Types.LimiterHomotopy.UpperLimit then uMax
                              else if homotopyType == Types.LimiterHomotopy.LowerLimit then uMin
                              else 0);
            if strict then
              if homotopyType == Types.LimiterHomotopy.NoHomotopy then
                y = smooth(0, noEvent(if u > uMax then uMax else if u < uMin then uMin else u));
              else
                y = homotopy(actual = smooth(0, noEvent(if u > uMax then uMax else if u < uMin then uMin else u)),
                             simplified=simplifiedExpr);
              end if;
            else
              if homotopyType == Types.LimiterHomotopy.NoHomotopy then
                y = smooth(0,if u > uMax then uMax else if u < uMin then uMin else u);
              else
                y = homotopy(actual = smooth(0,if u > uMax then uMax else if u < uMin then uMin else u),
                             simplified=simplifiedExpr);
              end if;
            end if;
            annotation (
              Documentation(info="<html>
<p>
The Limiter block passes its input signal as output signal
as long as the input is within the specified upper and lower
limits. If this is not the case, the corresponding limits are passed
as output.
</p>
<p>
The parameter <code>homotopyType</code> in the Advanced tab specifies the
simplified behaviour if homotopy-based initialization is used:
</p>
<ul>
<li><code>NoHomotopy</code>: the actual expression with limits is used</li>
<li><code>Linear</code>: a linear behaviour y = u is assumed (default option)</li>
<li><code>UpperLimit</code>: it is assumed that the output is stuck at the upper limit u = uMax</li>
<li><code>LowerLimit</code>: it is assumed that the output is stuck at the lower limit u = uMin</li>
</ul>
<p>
If it is known a priori in which region the input signal will be located, this option can help
a lot by removing one strong nonlinearity from the initialization problem.
</p>
</html>"),     Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{0,-90},{0,68}}, color={192,192,192}),
              Polygon(
                points={{0,90},{-8,68},{8,68},{0,90}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-90,0},{68,0}}, color={192,192,192}),
              Polygon(
                points={{90,0},{68,-8},{68,8},{90,0}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-80,-70},{-50,-70},{50,70},{80,70}}),
              Text(
                extent={{-150,-150},{150,-110}},
                textString="uMax=%uMax"),
              Line(
                visible=strict,
                points={{50,70},{80,70}},
                color={255,0,0}),
              Line(
                visible=strict,
                points={{-80,-70},{-50,-70}},
                color={255,0,0})}));
          end Limiter;
          annotation (
            Documentation(info="<html>
<p>
This package contains <strong>discontinuous</strong> and
<strong>non-differentiable, algebraic</strong> input/output blocks.
</p>
</html>",     revisions="<html>
<ul>
<li><em>October 21, 2002</em>
       by Christian Schweiger:<br>
       New block VariableLimiter added.</li>
<li><em>August 22, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized, based on an existing Dymola library
       of Dieter Moormann and Hilding Elmqvist.
</li>
</ul>
</html>"),     Icon(graphics={Line(points={{-80,-66},{-26,-66},{28,52},{88,52}},
                color={95,95,95})}));
    end Nonlinear;

    package Sources "Library of signal source blocks generating Real, Integer and Boolean signals"
      import Modelica.Blocks.Interfaces;
      extends Modelica.Icons.SourcesPackage;

      block RealExpression "Set output signal to a time varying Real expression"

        Modelica.Blocks.Interfaces.RealOutput y=0.0 "Value of Real output"
          annotation (Dialog(group="Time varying output signal"), Placement(
              transformation(extent={{100,-10},{120,10}})));

        annotation (Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}), graphics={
              Rectangle(
                extent={{-100,40},{100,-40}},
                fillColor={235,235,235},
                fillPattern=FillPattern.Solid,
                borderPattern=BorderPattern.Raised),
              Text(
                extent={{-96,15},{96,-15}},
                textString="%y"),
              Text(
                extent={{-150,90},{150,50}},
                textString="%name",
                textColor={0,0,255})}), Documentation(info="<html>
<p>
The (time varying) Real output signal of this block can be defined in its
parameter menu via variable <strong>y</strong>. The purpose is to support the
easy definition of Real expressions in a block diagram. For example,
in the y-menu the definition \"if time &lt; 1 then 0 else 1\" can be given in order
to define that the output signal is one, if time &ge; 1 and otherwise
it is zero. Note, that \"time\" is a built-in variable that is always
accessible and represents the \"model time\" and that
variable <strong>y</strong> is both a variable and a connector.
</p>
</html>"));

      end RealExpression;

      block Constant "Generate constant signal of type Real"
        parameter Real k(start=1) "Constant output value"
        annotation(Dialog(groupImage="modelica://Modelica/Resources/Images/Blocks/Sources/Constant.png"));
        extends Interfaces.SO;

      equation
        y = k;
        annotation (
          defaultComponentName="const",
          Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{-80,68},{-80,-80}}, color={192,192,192}),
              Polygon(
                points={{-80,90},{-88,68},{-72,68},{-80,90}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-90,-70},{82,-70}}, color={192,192,192}),
              Polygon(
                points={{90,-70},{68,-62},{68,-78},{90,-70}},
                lineColor={192,192,192},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid),
              Line(points={{-80,0},{80,0}}),
              Text(
                extent={{-150,-150},{150,-110}},
                textString="k=%k")}),
          Documentation(info="<html>
<p>
The Real output y is a constant signal:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Sources/Constant.png\"
     alt=\"Constant.png\">
</p>
</html>"));
      end Constant;

      block CombiTimeTable
        "Table look-up with respect to time and linear/periodic extrapolation methods (data from matrix/file)"
        import Modelica.Blocks.Tables.Internal;
        extends Modelica.Blocks.Interfaces.MO(final nout=max([size(columns, 1); size(offset, 1)]));
        parameter Boolean tableOnFile=false
          "= true, if table is defined on file or in function usertab"
          annotation (Dialog(group="Table data definition"));
        parameter Real table[:, :] = fill(0.0, 0, 2)
          "Table matrix (time = first column; e.g., table=[0, 0; 1, 1; 2, 4])"
          annotation (Dialog(group="Table data definition",enable=not tableOnFile));
        parameter String tableName="NoName"
          "Table name on file or in function usertab (see docu)"
          annotation (Dialog(group="Table data definition",enable=tableOnFile));
        parameter String fileName="NoName" "File where matrix is stored"
          annotation (Dialog(
            group="Table data definition",
            enable=tableOnFile,
            loadSelector(filter="Text files (*.txt);;MATLAB MAT-files (*.mat)",
                caption="Open file in which table is present")));
        parameter Boolean verboseRead=true
          "= true, if info message that file is loading is to be printed"
          annotation (Dialog(group="Table data definition",enable=tableOnFile));
        parameter Integer columns[:]=2:size(table, 2)
          "Columns of table to be interpolated"
          annotation (Dialog(group="Table data interpretation",
          groupImage="modelica://Modelica/Resources/Images/Blocks/Sources/CombiTimeTable.png"));
        parameter Modelica.Blocks.Types.Smoothness smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments
          "Smoothness of table interpolation"
          annotation (Dialog(group="Table data interpretation"));
        parameter Modelica.Blocks.Types.Extrapolation extrapolation=Modelica.Blocks.Types.Extrapolation.LastTwoPoints
          "Extrapolation of data outside the definition range"
          annotation (Dialog(group="Table data interpretation"));
        parameter SI.Time timeScale(
          min=Modelica.Constants.eps)=1 "Time scale of first table column"
          annotation (Dialog(group="Table data interpretation"), Evaluate=true);
        parameter Real offset[:]={0} "Offsets of output signals"
          annotation (Dialog(group="Table data interpretation"));
        parameter SI.Time startTime=0
          "Output = offset for time < startTime"
          annotation (Dialog(group="Table data interpretation"));
        parameter SI.Time shiftTime=startTime
          "Shift time of first table column"
          annotation (Dialog(group="Table data interpretation"));
        parameter Modelica.Blocks.Types.TimeEvents timeEvents=Modelica.Blocks.Types.TimeEvents.Always
          "Time event handling of table interpolation"
          annotation (Dialog(group="Table data interpretation", enable=smoothness == Modelica.Blocks.Types.Smoothness.LinearSegments));
        parameter Boolean verboseExtrapolation=false
          "= true, if warning messages are to be printed if time is outside the table definition range"
          annotation (Dialog(group="Table data interpretation", enable=extrapolation == Modelica.Blocks.Types.Extrapolation.LastTwoPoints or extrapolation == Modelica.Blocks.Types.Extrapolation.HoldLastPoint));
        final parameter SI.Time t_min=t_minScaled*timeScale
          "Minimum abscissa value defined in table";
        final parameter SI.Time t_max=t_maxScaled*timeScale
          "Maximum abscissa value defined in table";
        final parameter Real t_minScaled=Internal.getTimeTableTmin(tableID)
          "Minimum (scaled) abscissa value defined in table";
        final parameter Real t_maxScaled=Internal.getTimeTableTmax(tableID)
          "Maximum (scaled) abscissa value defined in table";
    protected
        final parameter Real p_offset[nout]=(if size(offset, 1) == 1 then ones(nout)*offset[1] else offset)
          "Offsets of output signals";
        parameter Modelica.Blocks.Types.ExternalCombiTimeTable tableID=
            Modelica.Blocks.Types.ExternalCombiTimeTable(
              if tableOnFile then tableName else "NoName",
              if tableOnFile and fileName <> "NoName" and not Modelica.Utilities.Strings.isEmpty(fileName) then fileName else "NoName",
              table,
              startTime/timeScale,
              columns,
              smoothness,
              extrapolation,
              shiftTime/timeScale,
              if smoothness == Modelica.Blocks.Types.Smoothness.LinearSegments then timeEvents elseif smoothness == Modelica.Blocks.Types.Smoothness.ConstantSegments then Modelica.Blocks.Types.TimeEvents.Always else Modelica.Blocks.Types.TimeEvents.NoTimeEvents,
              if tableOnFile then verboseRead else false) "External table object";
        discrete SI.Time nextTimeEvent(start=0, fixed=true)
          "Next time event instant";
        discrete Real nextTimeEventScaled(start=0, fixed=true)
          "Next scaled time event instant";
        Real timeScaled "Scaled time";
      equation
        if tableOnFile then
          assert(tableName <> "NoName",
            "tableOnFile = true and no table name given");
        else
          assert(size(table, 1) > 0 and size(table, 2) > 0,
            "tableOnFile = false and parameter table is an empty matrix");
        end if;

        if verboseExtrapolation and (
          extrapolation == Modelica.Blocks.Types.Extrapolation.LastTwoPoints or
          extrapolation == Modelica.Blocks.Types.Extrapolation.HoldLastPoint) then
          assert(noEvent(time >= t_min), "
Extrapolation warning: Time (="     + String(time) + ") must be greater or equal
than the minimum abscissa value t_min (="     + String(t_min) + ") defined in the table.
",     level=AssertionLevel.warning);
          assert(noEvent(time <= t_max), "
Extrapolation warning: Time (="     + String(time) + ") must be less or equal
than the maximum abscissa value t_max (="     + String(t_max) + ") defined in the table.
",     level=AssertionLevel.warning);
        end if;

        timeScaled = time/timeScale;
        when {time >= pre(nextTimeEvent), initial()} then
          nextTimeEventScaled = Internal.getNextTimeEvent(tableID, timeScaled);
          nextTimeEvent = if nextTimeEventScaled < Modelica.Constants.inf then nextTimeEventScaled*timeScale else Modelica.Constants.inf;
        end when;
        if smoothness == Modelica.Blocks.Types.Smoothness.ConstantSegments then
          for i in 1:nout loop
            y[i] = p_offset[i] + Internal.getTimeTableValueNoDer(tableID, i, timeScaled, nextTimeEventScaled, pre(nextTimeEventScaled));
          end for;
        elseif smoothness == Modelica.Blocks.Types.Smoothness.LinearSegments then
          for i in 1:nout loop
            y[i] = p_offset[i] + Internal.getTimeTableValueNoDer2(tableID, i, timeScaled, nextTimeEventScaled, pre(nextTimeEventScaled));
          end for;
        else
          for i in 1:nout loop
            y[i] = p_offset[i] + Internal.getTimeTableValue(tableID, i, timeScaled, nextTimeEventScaled, pre(nextTimeEventScaled));
          end for;
        end if;
        annotation (
          Documentation(info="<html>
<p>
This block generates an output signal y[:] by <strong>constant</strong>,
<strong>linear</strong> or <strong>cubic Hermite spline interpolation</strong>
in a table. The time points and function values are stored in a matrix
<strong>table[i,j]</strong>, where the first column table[:,1] contains the
time points and the other columns contain the data to be interpolated.
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Sources/CombiTimeTable.png\"
     alt=\"CombiTimeTable.png\">
</p>

<p>
Via parameter <strong>columns</strong> it can be defined which columns of the
table are interpolated. If, e.g., columns={2,4}, it is assumed that
2 output signals are present and that the first output is computed
by interpolation of column 2 and the second output is computed
by interpolation of column 4 of the table matrix.
The table interpolation has the following properties:
</p>
<ul>
<li>The interpolation interval is found by a binary search where the interval used in the
    last call is used as start interval.</li>
<li>The time points need to be <strong>strictly increasing</strong> for cubic Hermite
    spline interpolation, otherwise <strong>monotonically increasing</strong>.</li>
<li><strong>Discontinuities</strong> are allowed for (constant or) linear interpolation,
    by providing the same time point twice in the table.</li>
<li>Via parameter <strong>smoothness</strong> it is defined how the data is interpolated:
<blockquote><pre>
smoothness = 1: Linear interpolation
           = 2: Akima interpolation: Smooth interpolation by cubic Hermite
                splines such that der(y) is continuous, also if extrapolated.
           = 3: Constant segments
           = 4: Fritsch-Butland interpolation: Smooth interpolation by cubic
                Hermite splines such that y preserves the monotonicity and
                der(y) is continuous, also if extrapolated.
           = 5: Steffen interpolation: Smooth interpolation by cubic Hermite
                splines such that y preserves the monotonicity and der(y)
                is continuous, also if extrapolated.
           = 6: Modified Akima interpolation: Smooth interpolation by cubic
                Hermite splines such that der(y) is continuous, also if
                extrapolated. Additionally, overshoots and edge cases of the
                original Akima interpolation method are avoided.
</pre></blockquote></li>
<li>First and second <strong>derivatives</strong> are provided, with exception of the following two smoothness options.
<ol>
<li>No derivatives are provided for interpolation by constant segments.</li>
<li>No second derivative is provided for linear interpolation.<br>There is a design inconsistency, that it is possible
to model a signal consisting of constant segments using linear interpolation and duplicated sample points.
In contrast to interpolation by constant segments, the first derivative is provided as zero.</li>
</ol></li>
<li>Values <strong>outside</strong> of the table range, are computed by
    extrapolation according to the setting of parameter <strong>extrapolation</strong>:
<blockquote><pre>
extrapolation = 1: Hold the first or last value of the table,
                   if outside of the table scope.
              = 2: Extrapolate by using the derivative at the first/last table
                   points if outside of the table scope.
                   (If smoothness is LinearSegments or ConstantSegments
                   this means to extrapolate linearly through the first/last
                   two table points.).
              = 3: Periodically repeat the table data (periodical function).
              = 4: No extrapolation, i.e. extrapolation triggers an error
</pre></blockquote></li>
<li>If the table has only <strong>one row</strong>, no interpolation is performed and
    the table values of this row are just returned.</li>
<li>Via parameters <strong>shiftTime</strong> and <strong>offset</strong> the curve defined
    by the table can be shifted both in time and in the ordinate value.
    The time instants stored in the table are therefore <strong>relative</strong>
    to <strong>shiftTime</strong>.</li>
<li>If time &lt; startTime, no interpolation is performed and the offset
    is used as ordinate value for all outputs.</li>
<li>The table is implemented in a numerically sound way by
    generating <strong>time events</strong> at interval boundaries, in case of
    interpolation by linear segments.
    This generates continuously differentiable values for the integrator.
    Via parameter <strong>timeEvents</strong> it is defined how the time events are generated:
<blockquote><pre>
timeEvents = 1: Always generate time events at interval boundaries
           = 2: Generate time events at discontinuities (defined by duplicated sample points)
           = 3: No time events at interval boundaries
</pre></blockquote>
    For interpolation by constant segments time events are always generated at interval boundaries.
    For smooth interpolation by cubic Hermite splines no time events are generated at interval boundaries.</li>
<li>Via parameter <strong>timeScale</strong> the first column of the table array can
    be scaled, e.g., if the table array is given in hours (instead of seconds)
    <strong>timeScale</strong> shall be set to 3600.</li>
<li>For special applications it is sometimes needed to know the minimum
    and maximum time instant defined in the table as a parameter. For this
    reason parameters <strong>t_min</strong>/<strong>t_minScaled</strong> and
    <strong>t_max</strong>/<strong>t_maxScaled</strong> are provided and can be
    accessed from the outside of the table object. Whereas <strong>t_min</strong> and
    <strong>t_max</strong> define the scaled abscissa values (using parameter
    <strong>timeScale</strong>) in SI.Time, <strong>t_minScaled</strong> and
    <strong>t_maxScaled</strong> define the unitless original abscissa values of
    the table.</li>
</ul>
<p>
Example:
</p>
<blockquote><pre>
table = [0, 0;
         1, 0;
         1, 1;
         2, 4;
         3, 9;
         4, 16];
extrapolation = 2 (default), timeEvents = 2
If, e.g., time = 1.0, the output y =  0.0 (before event), 1.0 (after event)
    e.g., time = 1.5, the output y =  2.5,
    e.g., time = 2.0, the output y =  4.0,
    e.g., time = 5.0, the output y = 23.0 (i.e., extrapolation via last 2 points).
</pre></blockquote>
<p>
The table matrix can be defined in the following ways:
</p>
<ol>
<li>Explicitly supplied as <strong>parameter matrix</strong> \"table\",
    and the other parameters have the following values:
<blockquote><pre>
tableName is \"NoName\" or has only blanks,
fileName  is \"NoName\" or has only blanks.
</pre></blockquote></li>
<li><strong>Read</strong> from a <strong>file</strong> \"fileName\" where the matrix is stored as
    \"tableName\". Both text and MATLAB MAT-file format is possible.
    (The text format is described below).
    The MAT-file format comes in four different versions: v4, v6, v7 and v7.3.
    The library supports at least v4, v6 and v7 whereas v7.3 is optional.
    It is most convenient to generate the MAT-file from FreeMat or MATLAB&reg;
    by command
<blockquote><pre>
save tables.mat tab1 tab2 tab3
</pre></blockquote>
    or Scilab by command
<blockquote><pre>
savematfile tables.mat tab1 tab2 tab3
</pre></blockquote>
    when the three tables tab1, tab2, tab3 should be used from the model.<br>
    Note, a fileName can be defined as URI by using the helper function
    <a href=\"modelica://Modelica.Utilities.Files.loadResource\">loadResource</a>.</li>
<li>Statically stored in function \"usertab\" in file \"usertab.c\".
    The matrix is identified by \"tableName\". Parameter
    fileName = \"NoName\" or has only blanks. Row-wise storage is always to be
    preferred as otherwise the table is reallocated and transposed.</li>
</ol>
<p>
When the constant \"NO_FILE_SYSTEM\" is defined, all file I/O related parts of the
source code are removed by the C-preprocessor, such that no access to files takes place.
</p>
<p>
If tables are read from a text file, the file needs to have the
following structure (\"-----\" is not part of the file content):
</p>
<blockquote><pre>
-----------------------------------------------------
#1
double tab1(6,2)   # comment line
  0   0
  1   0
  1   1
  2   4
  3   9
  4  16
double tab2(6,2)   # another comment line
  0   0
  2   0
  2   2
  4   8
  6  18
  8  32
-----------------------------------------------------
</pre></blockquote>
<p>
Note, that the first two characters in the file need to be
\"#1\" (a line comment defining the version number of the file format).
Afterwards, the corresponding matrix has to be declared
with type (= \"double\" or \"float\"), name and actual dimensions.
Finally, in successive rows of the file, the elements of the matrix
have to be given. The elements have to be provided as a sequence of
numbers in row-wise order (therefore a matrix row can span several
lines in the file and need not start at the beginning of a line).
Numbers have to be given according to C syntax (such as 2.3, -2, +2.e4).
Number separators are spaces, tab (\\t), comma (,), or semicolon (;).
Several matrices may be defined one after another. Line comments start
with the hash symbol (#) and can appear everywhere.
Text files should either be ASCII or UTF-8 encoded, where UTF-8 encoded strings are only allowed in line comments and an optional UTF-8 BOM at the start of the text file is ignored.
Other characters, like trailing non comments, are not allowed in the file.
</p>
<p>
MATLAB is a registered trademark of The MathWorks, Inc.
</p>
</html>",     revisions="<html>
<p><strong>Release Notes:</strong></p>
<ul>
<li><em>April 09, 2013</em>
       by Thomas Beutlich:<br>
       Implemented as external object.</li>
<li><em>March 31, 2001</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Used CombiTableTime as a basis and added the
       arguments <strong>extrapolation, columns, startTime</strong>.
       This allows periodic function definitions.</li>
</ul>
</html>"),Icon(
          coordinateSystem(preserveAspectRatio=true,
            extent={{-100.0,-100.0},{100.0,100.0}}),
            graphics={
          Polygon(lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid,
            points={{-80.0,90.0},{-88.0,68.0},{-72.0,68.0},{-80.0,90.0}}),
          Line(points={{-80.0,68.0},{-80.0,-80.0}},
            color={192,192,192}),
          Line(points={{-90.0,-70.0},{82.0,-70.0}},
            color={192,192,192}),
          Polygon(lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid,
            points={{90.0,-70.0},{68.0,-62.0},{68.0,-78.0},{90.0,-70.0}}),
          Rectangle(lineColor={255,255,255},
            fillColor={255,215,136},
            fillPattern=FillPattern.Solid,
            extent={{-48.0,-50.0},{2.0,70.0}}),
          Line(points={{-48.0,-50.0},{-48.0,70.0},{52.0,70.0},{52.0,-50.0},{-48.0,-50.0},{-48.0,-20.0},{52.0,-20.0},{52.0,10.0},{-48.0,10.0},{-48.0,40.0},{52.0,40.0},{52.0,70.0},{2.0,70.0},{2.0,-51.0}})}));
      end CombiTimeTable;

      block BooleanTable
        "Generate a Boolean output signal based on a vector of time instants"

        parameter SI.Time table[:]={0,1}
          "Vector of time points. At every time point, the output y gets its opposite value (e.g., table={0,1})" annotation(Dialog(group="Table data definition"));
        parameter Boolean startValue=false
          "Start value of y. At time = table[1], y changes to 'not startValue'" annotation(Dialog(group="Table data interpretation",
          groupImage="modelica://Modelica/Resources/Images/Blocks/Sources/BooleanTable.png"));
        parameter Modelica.Blocks.Types.Extrapolation extrapolation=Modelica.Blocks.Types.Extrapolation.HoldLastPoint
          "Extrapolation of data outside the definition range" annotation(Dialog(group="Table data interpretation"));
        parameter SI.Time startTime=-Modelica.Constants.inf
          "Output = false for time < startTime" annotation(Dialog(group="Table data interpretation"));
        parameter SI.Time shiftTime=0
          "Shift time of table" annotation(Dialog(group="Table data interpretation"));

        extends Interfaces.partialBooleanSO;

        CombiTimeTable combiTimeTable(
          final table=if n > 0 then if startValue then [table[1], 1.0; table, {mod(i + 1, 2.0) for i in 1:n}] else [table[1], 0.0; table, {mod(i, 2.0) for i in 1:n}] else if startValue then [0.0, 1.0] else [0.0, 0.0],
          final smoothness=Modelica.Blocks.Types.Smoothness.ConstantSegments,
          final columns={2},
          final extrapolation=extrapolation,
          final startTime=startTime,
          final shiftTime=shiftTime) annotation(Placement(transformation(extent={{-30,-10},{-10,10}})));
        Modelica.Blocks.Math.RealToBoolean realToBoolean annotation(Placement(transformation(extent={{10,-10},{30,10}})));

    protected
          function isValidTable "Check if table is valid"
            extends Modelica.Icons.Function;
            input Real table[:] "Vector of time instants";
      protected
            Integer n=size(table, 1) "Number of table points";
          algorithm
            if n > 0 then
              // Check whether time values are strict monotonically increasing
              for i in 2:n loop
                assert(table[i] > table[i-1],
                  "Time values of table not strict monotonically increasing: table["
                   + String(i - 1) + "] = " + String(table[i - 1]) + ", table[" +
                  String(i) + "] = " + String(table[i]));
              end for;
            end if;
          end isValidTable;

          parameter Integer n=size(table, 1) "Number of table points";
      initial algorithm
          isValidTable(table);
      equation
          assert(extrapolation <> Modelica.Blocks.Types.Extrapolation.LastTwoPoints, "Unsuitable extrapolation setting.");
          connect(combiTimeTable.y[1], realToBoolean.u) annotation(Line(points={{-9,0},{8,0}}, color={0,0,127}));
          connect(realToBoolean.y, y) annotation(Line(points={{31,0},{110,0},{110,0}}, color={255,127,0}));
        annotation (
          Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                  100}}), graphics={Polygon(
                points={{-80,88},{-88,66},{-72,66},{-80,88}},
                lineColor={255,0,255},
                fillColor={255,0,255},
                fillPattern=FillPattern.Solid),
              Line(points={{-80,66},{-80,-82}}, color={255,0,255}),
              Line(points={{-90,-70},{72,-70}}, color={255,0,255}),
              Polygon(
                points={{90,-70},{68,-62},{68,-78},{90,-70}},
                lineColor={255,0,255},
                fillColor={255,0,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-18,70},{32,-50}},
                lineColor={255,255,255},
                fillColor={192,192,192},
                fillPattern=FillPattern.Solid), Line(points={{-18,-50},{-18,70},{32,
                    70},{32,-50},{-18,-50},{-18,-20},{32,-20},{32,10},{-18,10},{-18,
                    40},{32,40},{32,70},{32,70},{32,-51}})}),
          Documentation(info="<html>
<p>
The Boolean output y is a signal defined by parameter vector <strong>table</strong>.
In the vector time points are stored.
The table interpolation has the following properties:
</p>

<ul>
<li>At every time point, the output y
    changes its value to the negated value of the previous one.</li>
<li>Values <strong>outside</strong> of the table range, are computed by
    extrapolation according to the setting of parameter <strong>extrapolation</strong>:
<blockquote><pre>
extrapolation = 1: Hold the <strong>startValue</strong> or last value of the table,
                   if outside of the table scope.
              = 2: Extrapolate by using the derivative at the first/last table
                   points if outside of the table scope.
                   (This setting is not suitable and triggers an assert.)
              = 3: Periodically repeat the table data (periodical function).
              = 4: No extrapolation, i.e. extrapolation triggers an error
</pre></blockquote></li>
<li>Via parameter <strong>shiftTime</strong> the curve defined by the table can be shifted
    in time.
    The time instants stored in the table are therefore <strong>relative</strong>
    to <strong>shiftTime</strong>.</li>
<li>If time &lt; startTime, no interpolation is performed and <strong>false</strong>
    is used as ordinate value for the output.</li>
</ul>

<p>
<img src=\"modelica://Modelica/Resources/Images/Blocks/Sources/BooleanTable.png\"
     alt=\"BooleanTable.png\">
</p>

<p>
The precise semantics is:
</p>

<blockquote><pre>
<strong>if</strong> size(table,1) == 0 <strong>then</strong>
   y = startValue;
<strong>else</strong>
   //            time &lt; table[1]: y = startValue
   // table[1] &le; time &lt; table[2]: y = not startValue
   // table[2] &le; time &lt; table[3]: y = startValue
   // table[3] &le; time &lt; table[4]: y = not startValue
   // ...
<strong>end if</strong>;
</pre></blockquote>
</html>"));
      end BooleanTable;
      annotation (Documentation(info="<html>
<p>
This package contains <strong>source</strong> components, i.e., blocks which
have only output signals. These blocks are used as signal generators
for Real, Integer and Boolean signals.
</p>

<p>
All Real source signals (with the exception of the Constant source)
have at least the following two parameters:
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>offset</strong></td>
      <td>Value which is added to the signal</td>
  </tr>
  <tr><td><strong>startTime</strong></td>
      <td>Start time of signal. For time &lt; startTime,
                the output y is set to offset.</td>
  </tr>
</table>

<p>
The <strong>offset</strong> parameter is especially useful in order to shift
the corresponding source, such that at initial time the system
is stationary. To determine the corresponding value of offset,
usually requires a trimming calculation.
</p>
</html>",     revisions="<html>
<ul>
<li><em>October 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       Integer sources added. Step, TimeTable and BooleanStep slightly changed.</li>
<li><em>Nov. 8, 1999</em>
       by <a href=\"mailto:christoph@clauss-it.com\">Christoph Clau&szlig;</a>,
       <a href=\"mailto:Andre.Schneider@eas.iis.fraunhofer.de\">Andre.Schneider@eas.iis.fraunhofer.de</a>,
       <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       New sources: Exponentials, TimeTable. Trapezoid slightly enhanced
       (nperiod=-1 is an infinite number of periods).</li>
<li><em>Oct. 31, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       <a href=\"mailto:christoph@clauss-it.com\">Christoph Clau&szlig;</a>,
       <a href=\"mailto:Andre.Schneider@eas.iis.fraunhofer.de\">Andre.Schneider@eas.iis.fraunhofer.de</a>,
       All sources vectorized. New sources: ExpSine, Trapezoid,
       BooleanConstant, BooleanStep, BooleanPulse, SampleTrigger.
       Improved documentation, especially detailed description of
       signals in diagram layer.</li>
<li><em>June 29, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized a first version, based on an existing Dymola library
       of Dieter Moormann and Hilding Elmqvist.</li>
</ul>
</html>"));
    end Sources;

    package Tables "Library of blocks to interpolate in one and two-dimensional tables"
      extends Modelica.Icons.Package;

      package Internal "Internal external object definitions for table functions that should not be directly utilized by the user"
        extends Modelica.Icons.InternalPackage;

        pure function getTimeTableValue
          "Interpolate 1-dim. table where first column is time"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTimeTable tableID "External table object";
          input Integer icol "Column number";
          input Real timeIn "(Scaled) time value";
          discrete input Real nextTimeEvent "(Scaled) next time event in table";
          discrete input Real pre_nextTimeEvent "Pre-value of (scaled) next time event in table";
          output Real y "Interpolated value";
          external "C" y = ModelicaStandardTables_CombiTimeTable_getValue(tableID, icol, timeIn, nextTimeEvent, pre_nextTimeEvent)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
          annotation (derivative(
              noDerivative=nextTimeEvent,
              noDerivative=pre_nextTimeEvent) = getDerTimeTableValue);
        end getTimeTableValue;

        pure function getTimeTableValueNoDer
          "Interpolate 1-dim. table where first column is time (but do not provide a derivative function)"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTimeTable tableID "External table object";
          input Integer icol "Column number";
          input Real timeIn "(Scaled) time value";
          discrete input Real nextTimeEvent "(Scaled) next time event in table";
          discrete input Real pre_nextTimeEvent "Pre-value of (scaled) next time event in table";
          output Real y "Interpolated value";
          external "C" y = ModelicaStandardTables_CombiTimeTable_getValue(tableID, icol, timeIn, nextTimeEvent, pre_nextTimeEvent)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getTimeTableValueNoDer;

        pure function getTimeTableValueNoDer2
          "Interpolate 1-dim. table where first column is time (but do not provide a second derivative function)"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTimeTable tableID "External table object";
          input Integer icol "Column number";
          input Real timeIn "(Scaled) time value";
          discrete input Real nextTimeEvent "(Scaled) next time event in table";
          discrete input Real pre_nextTimeEvent "Pre-value of (scaled) next time event in table";
          output Real y "Interpolated value";
          external "C" y = ModelicaStandardTables_CombiTimeTable_getValue(tableID, icol, timeIn, nextTimeEvent, pre_nextTimeEvent)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
          annotation (derivative(
              noDerivative=nextTimeEvent,
              noDerivative=pre_nextTimeEvent) = getDerTimeTableValueNoDer);
        end getTimeTableValueNoDer2;

        pure function getDerTimeTableValue
          "Derivative of interpolated 1-dim. table where first column is time"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTimeTable tableID "External table object";
          input Integer icol "Column number";
          input Real timeIn "(Scaled) time value";
          discrete input Real nextTimeEvent "(Scaled) next time event in table";
          discrete input Real pre_nextTimeEvent "Pre-value of (scaled) next time event in table";
          input Real der_timeIn "Derivative of (scaled) time value";
          output Real der_y "Derivative of interpolated value";
          external "C" der_y = ModelicaStandardTables_CombiTimeTable_getDerValue(tableID, icol, timeIn, nextTimeEvent, pre_nextTimeEvent, der_timeIn)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
          annotation (derivative(
              order=2,
              noDerivative=nextTimeEvent,
              noDerivative=pre_nextTimeEvent) = getDer2TimeTableValue);
        end getDerTimeTableValue;

        pure function getDerTimeTableValueNoDer
          "Derivative of interpolated 1-dim. table where first column is time (but do not provide a derivative function)"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTimeTable tableID "External table object";
          input Integer icol "Column number";
          input Real timeIn "(Scaled) time value";
          discrete input Real nextTimeEvent "(Scaled) next time event in table";
          discrete input Real pre_nextTimeEvent "Pre-value of (scaled) next time event in table";
          input Real der_timeIn "Derivative of (scaled) time value";
          output Real der_y "Derivative of interpolated value";
          external "C" der_y = ModelicaStandardTables_CombiTimeTable_getDerValue(tableID, icol, timeIn, nextTimeEvent, pre_nextTimeEvent, der_timeIn)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getDerTimeTableValueNoDer;

        pure function getDer2TimeTableValue
          "Second derivative of interpolated 1-dim. table where first column is time"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTimeTable tableID "External table object";
          input Integer icol "Column number";
          input Real timeIn "(Scaled) time value";
          discrete input Real nextTimeEvent "(Scaled) next time event in table";
          discrete input Real pre_nextTimeEvent "Pre-value of (scaled) next time event in table";
          input Real der_timeIn "Derivative of (scaled) time value";
          input Real der2_timeIn "Second derivative of (scaled) time value";
          output Real der2_y "Second derivative of interpolated value";
          external "C" der2_y = ModelicaStandardTables_CombiTimeTable_getDer2Value(tableID, icol, timeIn, nextTimeEvent, pre_nextTimeEvent, der_timeIn, der2_timeIn)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getDer2TimeTableValue;

        pure function getTimeTableTmin
          "Return minimum abscissa value of 1-dim. table where first column is time"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTimeTable tableID "External table object";
          output Real timeMin "Minimum abscissa value in table";
          external "C" timeMin = ModelicaStandardTables_CombiTimeTable_minimumTime(tableID)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getTimeTableTmin;

        pure function getTimeTableTmax
          "Return maximum abscissa value of 1-dim. table where first column is time"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTimeTable tableID "External table object";
          output Real timeMax "Maximum abscissa value in table";
          external "C" timeMax = ModelicaStandardTables_CombiTimeTable_maximumTime(tableID)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getTimeTableTmax;

        pure function getNextTimeEvent
          "Return next time event value of 1-dim. table where first column is time"
          extends Modelica.Icons.Function;
          input Modelica.Blocks.Types.ExternalCombiTimeTable tableID "External table object";
          input Real timeIn "(Scaled) time value";
          output Real nextTimeEvent "(Scaled) next time event in table";
          external "C" nextTimeEvent = ModelicaStandardTables_CombiTimeTable_nextTimeEvent(tableID, timeIn)
            annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end getNextTimeEvent;
      end Internal;
      annotation (Documentation(info="<html>
<p>This package contains blocks for one- and two-dimensional interpolation in tables.</p>
<h4>Special interest topic: Statically stored tables for real-time simulation targets</h4>
<p>Especially for use on real-time platform targets (e.g., HIL-simulators) with <strong>no file system</strong>, it is possible to statically
store tables using a function &quot;usertab&quot; in a file conventionally named &quot;usertab.c&quot;. This can be more efficient than providing the tables as Modelica parameter arrays.</p>
<p>This is achieved by providing the tables in a specific structure as C-code and compiling that C-code together with the rest of the simulation model into a binary
that can be executed on the target platform. The &quot;Resources/Data/Tables/&quot; subdirectory of the MSL installation directory contains the files
<a href=\"modelica://Modelica/Resources/Data/Tables/usertab.c\">&quot;usertab.c&quot;</a> and <a href=\"modelica://Modelica/Resources/Data/Tables/usertab.h\">&quot;usertab.h&quot;</a>
that can be used as a template for own developments. While &quot;usertab.c&quot; would be typically used unmodified, the
&quot;usertab.h&quot; needs to adapted for the own needs.</p>
<p>In order to work it is necessary that the compiler pulls in the &quot;usertab.c&quot; file. Different Modelica tools might provide different mechanisms to do so.
Please consult the respective documentation/support for your Modelica tool.</p>
<p>A possible (though slightly makeshift) approach is to pull in the required files by utilizing a &quot;dummy&quot;-function that uses the Modelica external function
interface to include the required &quot;usertab.c&quot;. An example how this can be done is given below.</p>
<blockquote><pre>
model ExampleCTable \"Example utilizing the usertab.c interface\"
  extends Modelica.Icons.Example;
  parameter Real dummy(fixed=false) \"Dummy parameter\" annotation(HideResult=true);
  Modelica.Blocks.Tables.CombiTable1Dv table(tableOnFile=true, tableName=\"TestTable_1D_a\")
    annotation (Placement(transformation(extent={{-40,0},{-20,20}})));
  Modelica.Blocks.Sources.ContinuousClock clock
    annotation (Placement(transformation(extent={{-80,0},{-60,20}})));
protected
  encapsulated impure function getUsertab \"External dummy function to include \\\"usertab.c\\\"\"
    input Real dummy_u[:];
    output Real dummy_y;
    external \"C\" dummy_y = mydummyfunc(dummy_u);
    annotation(IncludeDirectory=\"modelica://Modelica/Resources/Data/Tables\",
           Include = \"#include \"usertab.c\"
double mydummyfunc(double* dummy_in) {
   return 0;
}
\");
  end getUsertab;
initial equation
  dummy = getUsertab(table.y);
equation
  connect(clock.y, table.u[1]) annotation (Line(points={{-59,10},{-42,10}}, color={0,0,127}));
  annotation (experiment(StartTime=0, StopTime=5), uses(Modelica(version=\"4.0.0\")));
end ExampleCTable;
</pre></blockquote>
</html>"),     Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
                {100,100}}), graphics={
            Rectangle(
              extent={{-76,-26},{80,-76}},
              lineColor={95,95,95},
              fillColor={235,235,235},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-76,24},{80,-26}},
              lineColor={95,95,95},
              fillColor={235,235,235},
              fillPattern=FillPattern.Solid),
            Rectangle(
              extent={{-76,74},{80,24}},
              lineColor={95,95,95},
              fillColor={235,235,235},
              fillPattern=FillPattern.Solid),
            Line(
              points={{-28,74},{-28,-76}},
              color={95,95,95}),
            Line(
              points={{24,74},{24,-76}},
              color={95,95,95})}));
    end Tables;

    package Types "Library of constants, external objects and types with choices, especially to build menus"
      extends Modelica.Icons.TypesPackage;

      type Smoothness = enumeration(
          LinearSegments "Linear interpolation of table points",
          ContinuousDerivative
            "Akima spline interpolation of table points (such that the first derivative is continuous)",
          ConstantSegments
            "Piecewise constant interpolation of table points (the value from the previous abscissa point is returned)",
          MonotoneContinuousDerivative1
            "Fritsch-Butland spline interpolation (such that the monotonicity is preserved and the first derivative is continuous)",
          MonotoneContinuousDerivative2
            "Steffen spline interpolation of table points (such that the monotonicity is preserved and the first derivative is continuous)",
          ModifiedContinuousDerivative
            "Modified Akima spline interpolation of table points (such that the first derivative is continuous and shortcomings of the original Akima method are avoided)")
        "Enumeration defining the smoothness of table interpolation";

        type Extrapolation = enumeration(
          HoldLastPoint
            "Hold the first/last table point outside of the table scope",
          LastTwoPoints
            "Extrapolate by using the derivative at the first/last table points outside of the table scope",
          Periodic "Repeat the table scope periodically",
          NoExtrapolation "Extrapolation triggers an error")
        "Enumeration defining the extrapolation of table interpolation";

        type TimeEvents = enumeration(
          Always "Always generate time events at interval boundaries",
          AtDiscontinuities "Generate time events at discontinuities (defined by duplicated sample points)",
          NoTimeEvents "No time events at interval boundaries")
        "Enumeration defining the time event handling of time table interpolation";

        type Init = enumeration(
          NoInit
            "No initialization (start values are used as guess values with fixed=false)",
          SteadyState
            "Steady state initialization (derivatives of states are zero)",
          InitialState "Initialization with initial states",
          InitialOutput
            "Initialization with initial outputs (and steady state of the states if possible)")
        "Enumeration defining initialization of a block" annotation (Evaluate=true,
        Documentation(info="<html>
  <p>The following initialization alternatives are available:</p>
  <dl>
    <dt><code><strong>NoInit</strong></code></dt>
      <dd>No initialization (start values are used as guess values with <code>fixed=false</code>)</dd>
    <dt><code><strong>SteadyState</strong></code></dt>
      <dd>Steady state initialization (derivatives of states are zero)</dd>
    <dt><code><strong>InitialState</strong></code></dt>
      <dd>Initialization with initial states</dd>
    <dt><code><strong>InitialOutput</strong></code></dt>
      <dd>Initialization with initial outputs (and steady state of the states if possible)</dd>
  </dl>
</html>"));

       type SimpleController = enumeration(
          P "P controller",
          PI "PI controller",
          PD "PD controller",
          PID "PID controller")
        "Enumeration defining P, PI, PD, or PID simple controller type" annotation (
         Evaluate=true);

      type LimiterHomotopy = enumeration(
          NoHomotopy "Homotopy is not used",
          Linear "Simplified model without limits",
          UpperLimit "Simplified model fixed at upper limit",
          LowerLimit "Simplified model fixed at lower limit")
        "Enumeration defining use of homotopy in limiter components" annotation (Evaluate=true);

      class ExternalCombiTimeTable
        "External object of 1-dim. table where first column is time"
        extends ExternalObject;

        function constructor "Initialize 1-dim. table where first column is time"
          extends Modelica.Icons.Function;
          input String tableName "Table name";
          input String fileName "File name";
          input Real table[:, :];
          input SI.Time startTime;
          input Integer columns[:];
          input Modelica.Blocks.Types.Smoothness smoothness;
          input Modelica.Blocks.Types.Extrapolation extrapolation;
          input SI.Time shiftTime=0.0;
          input Modelica.Blocks.Types.TimeEvents timeEvents=Modelica.Blocks.Types.TimeEvents.Always;
          input Boolean verboseRead=true "= true: Print info message; = false: No info message";
          output ExternalCombiTimeTable externalCombiTimeTable;
        external "C" externalCombiTimeTable = ModelicaStandardTables_CombiTimeTable_init2(
                fileName,
                tableName,
                table,
                size(table, 1),
                size(table, 2),
                startTime,
                columns,
                size(columns, 1),
                smoothness,
                extrapolation,
                shiftTime,
                timeEvents,
                verboseRead) annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end constructor;

        function destructor "Terminate 1-dim. table where first column is time"
          extends Modelica.Icons.Function;
          input ExternalCombiTimeTable externalCombiTimeTable;
        external "C" ModelicaStandardTables_CombiTimeTable_close(
            externalCombiTimeTable) annotation (IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStandardTables.h\"", Library={"ModelicaStandardTables", "ModelicaIO", "ModelicaMatIO", "zlib"});
        end destructor;

      end ExternalCombiTimeTable;
      annotation (Documentation(info="<html>
<p>
In this package <strong>types</strong>, <strong>constants</strong> and <strong>external objects</strong> are defined that are used
in library Modelica.Blocks. The types have additional annotation choices
definitions that define the menus to be built up in the graphical
user interface when the type is used as parameter in a declaration.
</p>
</html>"));
    end Types;

    package Icons "Icons for Blocks"
        extends Modelica.Icons.IconsPackage;

        partial block Block "Basic graphical layout of input/output block"

          annotation (
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                  100,100}}), graphics={Rectangle(
                extent={{-100,-100},{100,100}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{-150,150},{150,110}},
                textString="%name",
                textColor={0,0,255})}),
          Documentation(info="<html>
<p>
Block that has only the basic icon for an input/output
block (no declarations, no equations). Most blocks
of package Modelica.Blocks inherit directly or indirectly
from this block.
</p>
</html>"));

        end Block;

      partial block PartialBooleanBlock "Basic graphical layout of logical block"

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                fillColor={210,210,210},
                fillPattern=FillPattern.Solid,
                borderPattern=BorderPattern.Raised), Text(
                extent={{-150,150},{150,110}},
                textString="%name",
                textColor={0,0,255})}), Documentation(info="<html>
<p>
Block that has only the basic icon for an input/output,
Boolean block (no declarations, no equations) used especially
in the Blocks.Logical library.
</p>
</html>"));
      end PartialBooleanBlock;
    end Icons;
  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
        Rectangle(
          origin={0.0,35.1488},
          fillColor={255,255,255},
          extent={{-30.0,-20.1488},{30.0,20.1488}}),
        Rectangle(
          origin={0.0,-34.8512},
          fillColor={255,255,255},
          extent={{-30.0,-20.1488},{30.0,20.1488}}),
        Line(
          origin={-51.25,0.0},
          points={{21.25,-35.0},{-13.75,-35.0},{-13.75,35.0},{6.25,35.0}}),
        Polygon(
          origin={-40.0,35.0},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{10.0,0.0},{-5.0,5.0},{-5.0,-5.0}}),
        Line(
          origin={51.25,0.0},
          points={{-21.25,35.0},{13.75,35.0},{13.75,-35.0},{-6.25,-35.0}}),
        Polygon(
          origin={40.0,-35.0},
          pattern=LinePattern.None,
          fillPattern=FillPattern.Solid,
          points={{-10.0,0.0},{5.0,5.0},{5.0,-5.0}})}), Documentation(info="<html>
<p>
This library contains input/output blocks to build up block diagrams.
</p>

<dl>
<dt><strong>Main Author:</strong></dt>
<dd><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a><br>
    Deutsches Zentrum f&uuml;r Luft und Raumfahrt e. V. (DLR)<br>
    Oberpfaffenhofen<br>
    Postfach 1116<br>
    D-82230 Wessling<br>
    email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a><br></dd>
</dl>
<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>June 23, 2004</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Introduced new block connectors and adapted all blocks to the new connectors.
       Included subpackages Continuous, Discrete, Logical, Nonlinear from
       package ModelicaAdditions.Blocks.
       Included subpackage ModelicaAdditions.Table in Modelica.Blocks.Sources
       and in the new package Modelica.Blocks.Tables.
       Added new blocks to Blocks.Sources and Blocks.Logical.
       </li>
<li><em>October 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       New subpackage Examples, additional components.
       </li>
<li><em>June 20, 2000</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and
       Michael Tiller:<br>
       Introduced a replaceable signal type into
       Blocks.Interfaces.RealInput/RealOutput:
<blockquote><pre>
replaceable type SignalType = Real
</pre></blockquote>
       in order that the type of the signal of an input/output block
       can be changed to a physical type, for example:
<blockquote><pre>
Sine sin1(outPort(redeclare type SignalType=Modelica.Units.SI.Torque))
</pre></blockquote>
      </li>
<li><em>Sept. 18, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Renamed to Blocks. New subpackages Math, Nonlinear.
       Additional components in subpackages Interfaces, Continuous
       and Sources.</li>
<li><em>June 30, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized a first version, based on an existing Dymola library
       of Dieter Moormann and Hilding Elmqvist.</li>
</ul>
</html>"));
  end Blocks;

  package Fluid "Library of 1-dim. thermo-fluid flow models using the Modelica.Media media description"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;
    import Cv = Modelica.Units.Conversions;

    model System
      "System properties and default values (ambient, flow direction, initialization)"

      package Medium = Modelica.Media.Interfaces.PartialMedium
        "Medium model for default start values"
          annotation (choicesAllMatching = true);
      parameter SI.AbsolutePressure p_ambient=101325
        "Default ambient pressure"
        annotation(Dialog(group="Environment"));
      parameter SI.Temperature T_ambient=293.15
        "Default ambient temperature"
        annotation(Dialog(group="Environment"));
      parameter SI.Acceleration g=Modelica.Constants.g_n
        "Constant gravity acceleration"
        annotation(Dialog(group="Environment"));

      // Assumptions
      parameter Boolean allowFlowReversal = true
        "= false to restrict to design flow direction (port_a -> port_b)"
        annotation(Dialog(tab="Assumptions"), Evaluate=true);
      parameter Modelica.Fluid.Types.Dynamics energyDynamics=
        Modelica.Fluid.Types.Dynamics.DynamicFreeInitial
        "Default formulation of energy balances"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      parameter Modelica.Fluid.Types.Dynamics massDynamics=
        energyDynamics "Default formulation of mass balances"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      final parameter Modelica.Fluid.Types.Dynamics substanceDynamics=
        massDynamics "Default formulation of substance balances"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      final parameter Modelica.Fluid.Types.Dynamics traceDynamics=
        massDynamics "Default formulation of trace substance balances"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
      parameter Modelica.Fluid.Types.Dynamics momentumDynamics=
        Modelica.Fluid.Types.Dynamics.SteadyState
        "Default formulation of momentum balances, if options available"
        annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

      // Initialization
      parameter SI.MassFlowRate m_flow_start = 0
        "Default start value for mass flow rates"
        annotation(Dialog(tab = "Initialization"));
      parameter SI.AbsolutePressure p_start = p_ambient
        "Default start value for pressures"
        annotation(Dialog(tab = "Initialization"));
      parameter SI.Temperature T_start = T_ambient
        "Default start value for temperatures"
        annotation(Dialog(tab = "Initialization"));
      // Advanced
      parameter Boolean use_eps_Re = false
        "= true to determine turbulent region automatically using Reynolds number"
        annotation(Evaluate=true, Dialog(tab = "Advanced"));
      parameter SI.MassFlowRate m_flow_nominal = if use_eps_Re then 1 else 1e2*m_flow_small
        "Default nominal mass flow rate"
        annotation(Dialog(tab="Advanced", enable = use_eps_Re));
      parameter Real eps_m_flow(min=0) = 1e-4
        "Regularization of zero flow for |m_flow| < eps_m_flow*m_flow_nominal"
        annotation(Dialog(tab = "Advanced", enable = use_eps_Re));
      parameter SI.AbsolutePressure dp_small(min=0) = 1
        "Default small pressure drop for regularization of laminar and zero flow"
        annotation(Dialog(tab="Advanced", group="Classic", enable = not use_eps_Re));
      parameter SI.MassFlowRate m_flow_small(min=0) = 1e-2
        "Default small mass flow rate for regularization of laminar and zero flow"
        annotation(Dialog(tab = "Advanced", group="Classic", enable = not use_eps_Re));
    initial equation
      //assert(use_eps_Re, "*** Using classic system.m_flow_small and system.dp_small."
      //       + " They do not distinguish between laminar flow and regularization of zero flow."
      //       + " Absolute small values are error prone for models with local nominal values."
      //       + " Moreover dp_small can generally be obtained automatically."
      //       + " Please update the model to new system.use_eps_Re = true  (see system, Advanced tab). ***",
      //       level=AssertionLevel.warning);
      annotation (
        defaultComponentName="system",
        defaultComponentPrefixes="inner",
        missingInnerMessage="
Your model is using an outer \"system\" component but
an inner \"system\" component is not defined.
For simulation drag Modelica.Fluid.System into your model
to specify system properties.",
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                100}}), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              lineColor={0,0,255},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textColor={0,0,255},
              textString="%name"),
            Line(points={{-86,-30},{82,-30}}),
            Line(points={{-82,-68},{-52,-30}}),
            Line(points={{-48,-68},{-18,-30}}),
            Line(points={{-14,-68},{16,-30}}),
            Line(points={{22,-68},{52,-30}}),
            Line(points={{74,84},{74,14}}),
            Polygon(
              points={{60,14},{88,14},{74,-18},{60,14}},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{16,20},{60,-18}},
              textString="g"),
            Text(
              extent={{-90,82},{74,50}},
              textString="defaults"),
            Line(
              points={{-82,14},{-42,-20},{2,30}},
              thickness=0.5),
            Ellipse(
              extent={{-10,40},{12,18}},
              pattern=LinePattern.None,
              fillColor={255,0,0},
              fillPattern=FillPattern.Solid)}),
        Documentation(info="<html>
<p>
 A system component is needed in each fluid model to provide system-wide settings, such as ambient conditions and overall modeling assumptions.
 The system settings are propagated to the fluid models using the inner/outer mechanism.
</p>
<p>
 A model should never directly use system parameters.
 Instead a local parameter should be declared, which uses the global setting as default.
 The only exceptions are:</p>
 <ul>
  <li>the gravity system.g,</li>
  <li>the global system.eps_m_flow, which is used to define a local m_flow_small for the local m_flow_nominal:
      <blockquote><pre>m_flow_small = system.eps_m_flow*m_flow_nominal</pre></blockquote>
  </li>
 </ul>
<p>
 The global system.m_flow_small and system.dp_small are classic parameters.
 They do not distinguish between laminar flow and regularization of zero flow.
 Absolute small values are error prone for models with local nominal values.
 Moreover dp_small can generally be obtained automatically.
 Consider using the new system.use_eps_Re = true (see Advanced tab).
</p>
</html>"));
    end System;

    package Vessels "Devices for storing fluid"
        extends Modelica.Icons.VariantsPackage;

        model ClosedVolume
        "Volume of fixed size, closed to the ambient, with inlet/outlet ports"
        import Modelica.Constants.pi;

          // Mass and energy balance, ports
          extends Modelica.Fluid.Vessels.BaseClasses.PartialLumpedVessel(
            final fluidVolume = V,
            vesselArea = pi*(3/4*V)^(2/3),
            heatTransfer(surfaceAreas={4*pi*(3/4*V/pi)^(2/3)}));

          parameter SI.Volume V "Volume";

        equation
          Wb_flow = 0;
          for i in 1:nPorts loop
            vessel_ps_static[i] = medium.p;
          end for;

          annotation (defaultComponentName="volume",
            Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
                  100,100}}), graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillPattern=FillPattern.Sphere,
                fillColor={170,213,255}), Text(
                extent={{-150,12},{150,-18}},
                textString="V=%V")}),
          Documentation(info="<html>
<p>
Ideally mixed volume of constant size with two fluid ports and one medium model.
The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model if <code>use_portsData=false</code>.
Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected.
A spherical shape is assumed for the heat transfer area, with V=4/3*pi*r^3, A=4*pi*r^2.
Ideal heat transfer is assumed per default; the thermal port temperature is equal to the medium temperature.
</p>
<p>
If <code>use_portsData=true</code>, the port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe.
The hydraulic resistances <code>portsData.zeta_in</code> and <code>portsData.zeta_out</code> determine the dissipative pressure drop between volume and port depending on
the direction of mass flow. See <a href=\"modelica://Modelica.Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <em>[Idelchik, Handbook of Hydraulic Resistance, 2004]</em>.
</p>
</html>"));
        end ClosedVolume;

    model OpenTank "Simple tank with inlet/outlet ports"
        import Modelica.Constants.pi;

      // Tank properties
      SI.Height level(stateSelect=StateSelect.prefer, start=level_start_eps)
          "Level height of tank";
      SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";

      // Tank geometry
      parameter SI.Height height "Height of tank";
      parameter SI.Area crossArea "Area of tank";

      // Ambient
      parameter Medium.AbsolutePressure p_ambient=system.p_ambient
          "Tank surface pressure"
        annotation(Dialog(tab = "Assumptions", group = "Ambient"));
      parameter Medium.Temperature T_ambient=system.T_ambient
          "Tank surface Temperature"
        annotation(Dialog(tab = "Assumptions", group = "Ambient"));

      // Initialization
      parameter SI.Height level_start(min=0) = 0.5*height
          "Start value of tank level"
        annotation(Dialog(tab="Initialization"));

      // Mass and energy balance, ports
      extends Modelica.Fluid.Vessels.BaseClasses.PartialLumpedVessel(
        final fluidVolume = V,
        final fluidLevel = level,
        final fluidLevel_max = height,
        final vesselArea = crossArea,
        heatTransfer(surfaceAreas={crossArea+2*sqrt(crossArea*pi)*level}),
        final initialize_p = false,
        final p_start = p_ambient);

    protected
      final parameter SI.Height level_start_eps = max(level_start, Modelica.Constants.eps);

    equation
      // Total quantities
      V = crossArea*level "Volume of fluid";
      medium.p = p_ambient;

      // Source termsEnergy balance
      if Medium.singleState or energyDynamics == Types.Dynamics.SteadyState then
        Wb_flow = 0
            "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
      else
        Wb_flow = -p_ambient*der(V);
      end if;

      //Determine port properties
      for i in 1:nPorts loop
        vessel_ps_static[i] = max(0, level - portsData_height[i])*system.g*medium.d + p_ambient;
      end for;

    initial equation
      if massDynamics == Types.Dynamics.FixedInitial then
        level = level_start_eps;
      elseif massDynamics == Types.Dynamics.SteadyStateInitial then
        der(level) = 0;
      end if;

        annotation (defaultComponentName="tank",
          Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              initialScale=0.2), graphics={
              Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={255,255,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.VerticalCylinder),
              Rectangle(
                extent=DynamicSelect({{-100,-100},{100,10}}, {{-100,-100},{100,(-100
                     + 200*level/height)}}),
                fillColor={85,170,255},
                fillPattern=FillPattern.VerticalCylinder),
              Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}),
              Text(
                extent={{-95,60},{95,40}},
                textString="level ="),
              Text(
                extent={{-95,-24},{95,-44}},
                textString=DynamicSelect("%level_start", String(
                    level,
                    minimumLength=1,
                    significantDigits=2)))}),
          Documentation(info="<html>
<p>
Model of a tank that is open to the ambient at the fixed pressure
<code>p_ambient</code>.
</p>
<p>
The vector of connectors <strong>ports</strong> represents fluid ports at configurable heights, relative to the bottom of tank.
Fluid can flow either out of or in to each port.
</p>
The following assumptions are made:
<ul>
<li>The tank is filled with a single or multiple-substance medium having a density higher than the density of the ambient medium.</li>
<li>The fluid has uniform density, temperature and mass fractions</li>
<li>No liquid is leaving the tank through the open top; the simulation breaks with an assertion if the liquid level growths over the height.</li>
</ul>
<p>
The port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe.
The hydraulic resistances <code>portsData.zeta_in</code> and <code>portsData.zeta_out</code> determine the dissipative pressure drop between tank and port depending on
the direction of mass flow. See <a href=\"modelica://Modelica.Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <em>[Idelchik, Handbook of Hydraulic Resistance, 2004]</em>.
</p>
<p>
With the setting <code>use_portsData=false</code>, the port pressure represents the static head
at the height of the respective port.
The relationship between pressure drop and mass flow rate at the port must then be provided by connected components;
Heights of ports as well as kinetic and potential energy of fluid entering or leaving are not taken into account anymore.
</p>
</html>",     revisions="<html>
<ul>
<li><em>Dec. 12, 2008</em> by R&uuml;diger Franke: move port definitions
   to BaseClasses.PartialLumpedVessel; also use energy and mass balance from common base class</li>
<li><em>Dec. 8, 2008</em> by Michael Wetter (LBNL):<br>
Implemented trace substances.</li>
<li><em>Jan. 6, 2006</em> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund),
   Martin Otter (DLR):<br>
   Implementation based on former tank model.</li>
<li><em>Oct. 29, 2007</em> by Carsten Heinrich (ILK Dresden):<br>
Adapted to the new fluid library interfaces:
<ul> <li>FluidPorts_b is used instead of FluidPort_b (due to it is defined as an array of ports)</li>
    <li>Port name changed from port to ports</li></ul>Updated documentation.</li>
<li><em>Apr. 25, 2006</em> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
Limitation to bottom ports only, added inlet and outlet loss factors.</li>
</ul>
</html>"));
    end OpenTank;

      package BaseClasses "Base classes used in the Vessels package (only of interest to build new component models)"
        extends Modelica.Icons.BasesPackage;

          partial model PartialLumpedVessel
          "Lumped volume with a vector of fluid ports and replaceable heat transfer model"
            extends Modelica.Fluid.Interfaces.PartialLumpedVolume;

            // Port definitions
            parameter Integer nPorts=0 "Number of ports"
              annotation(Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));
            VesselFluidPorts_b ports[nPorts](redeclare each package Medium = Medium)
            "Fluid inlets and outlets"
              annotation (Placement(transformation(extent={{-40,-10},{40,10}},
                origin={0,-100})));

            // Port properties
            parameter Boolean use_portsData=true
            "= false to neglect pressure loss and kinetic energy"
              annotation(Evaluate=true, Dialog(tab="General",group="Ports"));
            parameter Modelica.Fluid.Vessels.BaseClasses.VesselPortsData[if use_portsData then nPorts else 0]
            portsData "Data of inlet/outlet ports"
              annotation(Dialog(tab="General",group="Ports",enable= use_portsData));

            parameter Medium.MassFlowRate m_flow_nominal = if system.use_eps_Re then system.m_flow_nominal else 1e2*system.m_flow_small
            "Nominal value for mass flow rates in ports"
              annotation(Dialog(tab="Advanced", group="Port properties"));
            parameter SI.MassFlowRate m_flow_small(min=0) = if system.use_eps_Re then system.eps_m_flow*m_flow_nominal else system.m_flow_small
            "Regularization range at zero mass flow rate"
              annotation(Dialog(tab="Advanced", group="Port properties"));
            parameter Boolean use_Re = system.use_eps_Re
            "= true, if turbulent region is defined by Re, otherwise by m_flow_small"
              annotation(Dialog(tab="Advanced", group="Port properties"), Evaluate=true);

            Medium.EnthalpyFlowRate ports_H_flow[nPorts];
            Medium.MassFlowRate ports_mXi_flow[nPorts,Medium.nXi];
            Medium.MassFlowRate[Medium.nXi] sum_ports_mXi_flow
            "Substance mass flows through ports";
            Medium.ExtraPropertyFlowRate ports_mC_flow[nPorts,Medium.nC];
            Medium.ExtraPropertyFlowRate[Medium.nC] sum_ports_mC_flow
            "Trace substance mass flows through ports";

            // Heat transfer through boundary
            parameter Boolean use_HeatTransfer = false
            "= true to use the HeatTransfer model"
                annotation (Dialog(tab="Assumptions", group="Heat transfer"));
            replaceable model HeatTransfer =
                Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer
              constrainedby Modelica.Fluid.Vessels.BaseClasses.HeatTransfer.PartialVesselHeatTransfer
            "Wall heat transfer"
                annotation (Dialog(tab="Assumptions", group="Heat transfer",enable=use_HeatTransfer),choicesAllMatching=true);
            HeatTransfer heatTransfer(
              redeclare final package Medium = Medium,
              final n=1,
              final states = {medium.state},
              final use_k = use_HeatTransfer)
                annotation (Placement(transformation(
                  extent={{-10,-10},{30,30}},
                  rotation=90,
                  origin={-50,-10})));
            Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort if use_HeatTransfer
              annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

            // Conservation of kinetic energy
            Medium.Density[nPorts] portInDensities
            "Densities of the fluid at the device boundary";
            SI.Velocity[nPorts] portVelocities
            "Velocities of fluid flow at device boundary";
            SI.EnergyFlowRate[nPorts] ports_E_flow
            "Flow of kinetic and potential energy at device boundary";

            // Note: should use fluidLevel_start - portsData.height
            Real[nPorts] s(each start = fluidLevel_max)
            "Curve parameters for port flows vs. port pressures; for further details see, Modelica Tutorial: Ideal switching devices";
            Real[nPorts] ports_penetration
            "Penetration of port with fluid, depending on fluid level and port diameter";

            // treatment of pressure losses at ports
            SI.Area[nPorts] portAreas = {Modelica.Constants.pi/4*portsData_diameter[i]^2 for i in 1:nPorts};
            Medium.AbsolutePressure[nPorts] vessel_ps_static
            "Static pressures inside the vessel at the height of the corresponding ports, zero flow velocity";

            // determination of turbulent region
            constant SI.ReynoldsNumber Re_turbulent = 100 "cf. suddenExpansion";
            SI.MassFlowRate[nPorts] m_flow_turbulent;

      protected
            input SI.Height fluidLevel = 0
            "Level of fluid in the vessel for treating heights of ports";
            parameter SI.Height fluidLevel_max = 1
            "Maximum level of fluid in the vessel";
            parameter SI.Area vesselArea = Modelica.Constants.inf
            "Area of the vessel used to relate to cross flow area of ports";

            // Treatment of use_portsData=false to neglect portsData and to not require its specification either in this case.
            // Remove portsData conditionally if use_portsData=false. Simplify their use in model equations by always
            // providing portsData_diameter and portsData_height, independent of the use_portsData setting.
            // Note: this moreover serves as work-around if a tool does not support a zero sized portsData record.
            Modelica.Blocks.Interfaces.RealInput[nPorts]
            portsData_diameter_internal = portsData.diameter if use_portsData and nPorts > 0;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height_internal = portsData.height if use_portsData and nPorts > 0;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in_internal = portsData.zeta_in if use_portsData and nPorts > 0;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out_internal = portsData.zeta_out if use_portsData and nPorts > 0;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_diameter;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in;
            Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out;
            Modelica.Blocks.Interfaces.BooleanInput[nPorts] regularFlow(each start=true);
            Modelica.Blocks.Interfaces.BooleanInput[nPorts] inFlow(each start=false);

          equation
            mb_flow = sum(ports.m_flow);
            mbXi_flow = sum_ports_mXi_flow;
            mbC_flow  = sum_ports_mC_flow;
            Hb_flow = sum(ports_H_flow) + sum(ports_E_flow);
            Qb_flow = heatTransfer.Q_flows[1];

            // Only one connection allowed to a port to avoid unwanted ideal mixing
            for i in 1:nPorts loop
              assert(cardinality(ports[i]) <= 1,"
each ports[i] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts to add an additional port.
");         end for;
            // Check for correct solution
            assert(fluidLevel <= fluidLevel_max, "Vessel is overflowing (fluidLevel > fluidLevel_max = " + String(fluidLevel) + ")");
            assert(fluidLevel > -1e-6*fluidLevel_max, "Fluid level (= " + String(fluidLevel) + ") is below zero meaning that the solution failed.");

            // Boundary conditions

            // treatment of conditional portsData
            connect(portsData_diameter, portsData_diameter_internal);
            connect(portsData_height, portsData_height_internal);
            connect(portsData_zeta_in, portsData_zeta_in_internal);
            connect(portsData_zeta_out, portsData_zeta_out_internal);
            if not use_portsData then
              portsData_diameter = zeros(nPorts);
              portsData_height = zeros(nPorts);
              portsData_zeta_in = zeros(nPorts);
              portsData_zeta_out = zeros(nPorts);
            end if;

            // actual definition of port variables
            for i in 1:nPorts loop
              portInDensities[i] = Medium.density(Medium.setState_phX(vessel_ps_static[i], inStream(ports[i].h_outflow), inStream(ports[i].Xi_outflow)));
              if use_portsData then
                // dp = 0.5*zeta*d*v*|v|
                // Note: assume vessel_ps_static for portVelocities to avoid algebraic loops for ports.p
                portVelocities[i] = smooth(0, ports[i].m_flow/portAreas[i]/Medium.density(Medium.setState_phX(vessel_ps_static[i], actualStream(ports[i].h_outflow), actualStream(ports[i].Xi_outflow))));
                // Note: the penetration should not go too close to zero as this would prevent a vessel from running empty
                ports_penetration[i] = Utilities.regStep(fluidLevel - portsData_height[i] - 0.1*portsData_diameter[i], 1, 1e-3, 0.1*portsData_diameter[i]);
                m_flow_turbulent[i]=if not use_Re then m_flow_small else
                  max(m_flow_small, (Modelica.Constants.pi/8)*portsData_diameter[i]
                                     *(Medium.dynamicViscosity(Medium.setState_phX(vessel_ps_static[i], inStream(ports[i].h_outflow), inStream(ports[i].Xi_outflow)))
                                       + Medium.dynamicViscosity(medium.state))*Re_turbulent);
              else
                // an infinite port diameter is assumed
                portVelocities[i] = 0;
                ports_penetration[i] = 1;
                m_flow_turbulent[i] = Modelica.Constants.inf;
              end if;

              // fluid flow through ports
              regularFlow[i] = fluidLevel >= portsData_height[i];
              inFlow[i]      = not regularFlow[i] and (s[i] > 0 or portsData_height[i] >= fluidLevel_max);
              if regularFlow[i] then
                // regular operation: fluidLevel is above ports[i]
                // Note: >= covers default values of zero as well
                if use_portsData then
                  /* Without regularization
                 ports[i].p = vessel_ps_static[i] + 0.5*ports[i].m_flow^2/portAreas[i]^2
                              * noEvent(if ports[i].m_flow>0 then zeta_in[i]/portInDensities[i] else -zeta_out[i]/medium.d);
              */

                  ports[i].p = vessel_ps_static[i] + (0.5/portAreas[i]^2*Utilities.regSquare2(ports[i].m_flow, m_flow_turbulent[i],
                                    (portsData_zeta_in[i] - 1 + portAreas[i]^2/vesselArea^2)/portInDensities[i]*ports_penetration[i],
                                    (portsData_zeta_out[i] + 1 - portAreas[i]^2/vesselArea^2)/medium.d/ports_penetration[i]));
                  /*
                // alternative formulation m_flow=f(dp); not allowing the ideal portsData_zeta_in[i]=1 though
                ports[i].m_flow = smooth(2, portAreas[i]*Utilities.regRoot2(ports[i].p - vessel_ps_static[i], dp_small,
                                       2*portInDensities[i]/portsData_zeta_in[i],
                                       2*medium.d/portsData_zeta_out[i]));
              */
                else
                  ports[i].p = vessel_ps_static[i];
                end if;
                s[i] = fluidLevel - portsData_height[i];

              elseif inFlow[i] then
                // ports[i] is above fluidLevel and has inflow
                ports[i].p = vessel_ps_static[i];
                s[i] = ports[i].m_flow;

              else
                // ports[i] is above fluidLevel, preventing outflow
                ports[i].m_flow = 0;
                s[i] = (ports[i].p - vessel_ps_static[i])/Medium.p_default*(portsData_height[i] - fluidLevel);
              end if;

              ports[i].h_outflow  = medium.h;
              ports[i].Xi_outflow = medium.Xi;
              ports[i].C_outflow  = C;

              ports_H_flow[i] = ports[i].m_flow * actualStream(ports[i].h_outflow)
              "Enthalpy flow";
              ports_E_flow[i] = ports[i].m_flow*(0.5*portVelocities[i]*portVelocities[i] + system.g*portsData_height[i])
              "Flow of kinetic and potential energy";
              ports_mXi_flow[i,:] = ports[i].m_flow * actualStream(ports[i].Xi_outflow)
              "Component mass flow";
              ports_mC_flow[i,:]  = ports[i].m_flow * actualStream(ports[i].C_outflow)
              "Trace substance mass flow";
            end for;

            for i in 1:Medium.nXi loop
              sum_ports_mXi_flow[i] = sum(ports_mXi_flow[:,i]);
            end for;

            for i in 1:Medium.nC loop
              sum_ports_mC_flow[i]  = sum(ports_mC_flow[:,i]);
            end for;

            connect(heatPort, heatTransfer.heatPorts[1]) annotation (Line(
                points={{-100,0},{-87,0},{-87,0},{-74,0}}, color={191,0,0}));
           annotation (
            Documentation(info="<html>
<p>
This base class extends PartialLumpedVolume with a vector of fluid ports and a replaceable wall HeatTransfer model.
</p>
<p>
The following modeling assumption are made:</p>
<ul>
<li>homogeneous medium, i.e., phase separation is not taken into account,</li>
<li>no kinetic energy in the fluid, i.e., kinetic energy dissipates into the internal energy,</li>
<li>pressure loss definitions at vessel ports assume incompressible fluid,</li>
<li>outflow of ambient media is prevented at each port assuming check valve behavior.
    If <code>fluidlevel &lt; portsData_height[i]</code> and <code>ports[i].p &lt; vessel_ps_static[i]</code> mass flow at the port is set to 0.</li>
</ul>
<p>
Each port has a (hydraulic) diameter and a height above the bottom of the vessel, which can be configured using the <strong><code>portsData</code></strong> record.
Alternatively the impact of port geometries can be neglected with <code>use_portsData=false</code>. This might be useful for early
design studies. Note that this means to assume an infinite port diameter at the bottom of the vessel.
Pressure drops and heights of the ports as well as kinetic and potential energy fluid entering or leaving the vessel are neglected then.
</p>
<p>
The following variables need to be defined by an extending model:
</p>
<ul>
<li><code>input fluidVolume</code>, the volume of the fluid in the vessel,</li>
<li><code>vessel_ps_static[nPorts]</code>, the static pressures inside the vessel at the height of the corresponding ports, at zero flow velocity, and</li>
<li><code>Wb_flow</code>, work term of the energy balance, e.g., p*der(V) if the volume is not constant or stirrer power.</li>
</ul>
<p>
An extending model should define:
</p>
<ul>
<li><code>parameter vesselArea</code> (default: Modelica.Constants.inf m2), the area of the vessel, to be related to cross flow areas of the ports for the consideration of dynamic pressure effects.</li>
</ul>
<p>
Optionally the fluid level may vary in the vessel, which effects the flow through the ports at configurable <code>portsData_height[nPorts]</code>.
This is why an extending model with varying fluid level needs to define:
</p>
<ul>
<li><code>input fluidLevel (default: 0m)</code>, the level the fluid in the vessel, and</li>
<li><code>parameter fluidLevel_max (default: 1m)</code>, the maximum level that must not be exceeded. Ports at or above fluidLevel_max can only receive inflow.</li>
</ul>
<p>
An extending model should not access the <code>portsData</code> record defined in the configuration dialog,
as an access to <code>portsData</code> may fail for <code>use_portsData=false</code> or <code>nPorts=0</code>.
</p>
<p>
Instead the predefined variables
</p>
<ul>
<li><code>portsData_diameter[nPorts]</code>,</li>
<li><code>portsData_height[nPorts]</code>,</li>
<li><code>portsData_zeta_in[nPorts]</code>, and</li>
<li><code>portsData_zeta_out[nPorts]</code></li>
</ul>
<p>
should be used if these values are needed.
</p>
</html>",           revisions="<html>
<ul>
<li><em>Jan. 2009</em> by R&uuml;diger Franke: extended with
   <ul><li>portsData record and threat configurable port heights,</li>
       <li>consideration of kinetic and potential energy of fluid entering or leaving in energy balance</li>
   </ul>
</li>
<li><em>Dec. 2008</em> by R&uuml;diger Franke: derived from OpenTank, in order to make general use of configurable port diameters</li>
</ul>
</html>"),    Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                    {100,100}}), graphics={Text(
                  extent={{-150,110},{150,150}},
                  textString="%name",
                  textColor={0,0,255})}));
          end PartialLumpedVessel;

      package HeatTransfer "HeatTransfer models for vessels"
        extends Modelica.Icons.Package;

        partial model PartialVesselHeatTransfer
            "Base class for vessel heat transfer models"
          extends Modelica.Fluid.Interfaces.PartialHeatTransfer;

          annotation(Documentation(info="<html>
Base class for vessel heat transfer models.
</html>"),    Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                      {100,100}}), graphics={Ellipse(
                    extent={{-60,64},{60,-56}},
                    fillPattern=FillPattern.Sphere,
                    fillColor={232,0,0}), Text(
                    extent={{-38,26},{40,-14}},
                    textString="%name")}));
        end PartialVesselHeatTransfer;

        model IdealHeatTransfer
            "IdealHeatTransfer: Ideal heat transfer without thermal resistance"
          extends PartialVesselHeatTransfer;

        equation
          Ts = heatPorts.T;

          annotation(Documentation(info="<html>
Ideal heat transfer without thermal resistance.
</html>"));
        end IdealHeatTransfer;
        annotation (Documentation(info="<html>
Heat transfer correlations for pipe models
</html>"));
      end HeatTransfer;

        record VesselPortsData "Data to describe inlet/outlet ports at vessels:
    diameter -- Inner (hydraulic) diameter of inlet/outlet port
    height -- Height over the bottom of the vessel
    zeta_out -- Hydraulic resistance out of vessel, default 0.5 for small diameter mounted flush with the wall
    zeta_in -- Hydraulic resistance into vessel, default 1.04 for small diameter mounted flush with the wall"
              extends Modelica.Icons.Record;
          parameter SI.Diameter diameter
            "Inner (hydraulic) diameter of inlet/outlet port";
          parameter SI.Height height = 0 "Height over the bottom of the vessel";
          parameter Real zeta_out(min=0)=0.5
            "Hydraulic resistance out of vessel, default 0.5 for small diameter mounted flush with the wall";
          parameter Real zeta_in(min=0)=1.04
            "Hydraulic resistance into vessel, default 1.04 for small diameter mounted flush with the wall";
          annotation (preferredView="info", Documentation(info="<html>
<h4>Vessel Port Data</h4>
<p>
This record describes the <strong>ports</strong> of a <strong>vessel</strong>. The variables in it are mostly self-explanatory (see list below); only the &zeta;
loss factors are discussed further. All data is quoted from Idelchik (1994).
</p>

<h4>Outlet Coefficients</h4>

<p>
If a <strong>straight pipe with constant cross section is mounted flush with the wall</strong>, its outlet pressure loss coefficient will be <code>&zeta; = 0.5</code> (Idelchik, p. 160, Diagram 3-1, paragraph 2).
</p>
<p>
If a <strong>straight pipe with constant cross section is mounted into a vessel such that the entrance into it is at a distance</strong> <code>b</code> from the wall (inside) the following table can be used. Herein, &delta; is the tube wall thickness (Idelchik, p. 160, Diagram 3-1, paragraph 1).
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <caption align=\"bottom\">Pressure loss coefficients for outlets, entrance at a distance from wall</caption>
  <tr>
    <td></td> <td>   </td><th colspan=\"5\" align=\"center\"> b / D_hyd  </th>
  </tr>
  <tr>
    <td></td> <td>   </td><th> 0.000 </th><th> 0.005 </th><th> 0.020 </th><th> 0.100 </th><th> 0.500-&#8734; </th>
  </tr>
  <tr>
     <th rowspan=\"5\" valign=\"middle\">&delta; / D_hyd</th> <th> 0.000 </th><td> 0.50 </td><td> 0.63  </td><td> 0.73  </td><td> 0.86  </td><td>      1.00     </td>
  </tr>
  <tr>
              <th> 0.008 </th><td> 0.50 </td><td> 0.55  </td><td> 0.62  </td><td> 0.74  </td><td>      0.88     </td>
  </tr>
  <tr>
              <th> 0.016 </th><td> 0.50 </td><td> 0.51  </td><td> 0.55  </td><td> 0.64  </td><td>      0.77     </td>
  </tr>
  <tr>
              <th> 0.024 </th><td> 0.50 </td><td> 0.50  </td><td> 0.52  </td><td> 0.58  </td><td>      0.68     </td>
  </tr>
  <tr>
              <th> 0.040 </th><td> 0.50 </td><td> 0.50  </td><td> 0.51  </td><td> 0.51  </td><td>      0.54     </td>
  </tr>
</table>

<p>
If a <strong>straight pipe with a circular bellmouth inlet (collector) without baffle is mounted flush with the wall</strong> then its pressure loss coefficient can be established from the following table. Herein, r is the radius of the bellmouth inlet surface (Idelchik, p. 164 f., Diagram 3-4, paragraph b)
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <caption align=\"bottom\">Pressure loss coefficients for outlets, bellmouth flush with wall</caption>
  <tr>
    <td></td> <th colspan=\"6\" align=\"center\"> r / D_hyd  </th>
  </tr>
  <tr>
    <td></td> <th> 0.01 </th><th> 0.03 </th><th> 0.05 </th><th> 0.08 </th><th> 0.16 </th><th>&ge;0.20</th>
  </tr>
  <tr>
     <th>&zeta;</th> <td> 0.44 </td><td> 0.31 </td><td> 0.22  </td><td> 0.15  </td><td> 0.06  </td><td>      0.03     </td>
  </tr>
</table>

<p>
If a <strong>straight pipe with a circular bellmouth inlet (collector) without baffle is mounted at a distance from a wall</strong> then its pressure loss coefficient can be established from the following table. Herein, r is the radius of the bellmouth inlet surface (Idelchik, p. 164 f., Diagram 3-4, paragraph a)
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <caption align=\"bottom\">Pressure loss coefficients for outlets, bellmouth at a distance of wall</caption>
  <tr>
    <td></td> <th colspan=\"6\" align=\"center\"> r / D_hyd  </th>
  </tr>
  <tr>
    <td></td> <th> 0.01 </th><th> 0.03 </th><th> 0.05 </th><th> 0.08 </th><th> 0.16 </th><th>&ge;0.20</th>
  </tr>
  <tr>
     <th>&zeta;</th> <td> 0.87 </td><td> 0.61 </td><td> 0.40  </td><td> 0.20  </td><td> 0.06  </td><td>      0.03     </td>
  </tr>
</table>

<h4>Inlet Coefficients</h4>

<p>
If a <strong>straight pipe with constant circular cross section is mounted flush with the wall</strong>, its vessel inlet pressure loss coefficient will be according to the following table (Idelchik, p. 209 f., Diagram 4-2 with <code>A_port/A_vessel = 0</code> and Idelchik, p. 640, Diagram 11-1, graph a). According to the text, <code>m = 9</code> is appropriate for fully developed turbulent flow.
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <caption align=\"bottom\">Pressure loss coefficients for inlets, circular tube flush with wall</caption>
  <tr>
    <td></td> <th colspan=\"6\" align=\"center\"> m  </th>
  </tr>
  <tr>
    <td></td> <th> 1.0 </th><th> 2.0 </th><th> 3.0 </th><th> 4.0 </th><th> 7.0 </th><th>9.0</th>
  </tr>
  <tr>
     <th>&zeta;</th> <td> 2.70 </td><td> 1.50 </td><td> 1.25  </td><td> 1.15  </td><td> 1.06  </td><td>      1.04     </td>
  </tr>
</table>

<p>
For larger port diameters, relative to the area of the vessel, the inlet pressure loss coefficient will be according to the following table (Idelchik, p. 209 f., Diagram 4-2 with <code>m = 7</code>).
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <caption align=\"bottom\">Pressure loss coefficients for inlets, circular tube flush with wall</caption>
  <tr>
    <td></td> <th colspan=\"6\" align=\"center\"> A_port / A_vessel  </th>
  </tr>
  <tr>
    <td></td> <th> 0.0 </th><th> 0.1 </th><th> 0.2 </th><th> 0.4 </th><th> 0.6 </th><th>0.8</th>
  </tr>
  <tr>
     <th>&zeta;</th> <td> 1.04 </td><td> 0.84 </td><td> 0.67  </td><td> 0.39  </td><td> 0.18  </td><td>      0.06     </td>
  </tr>
</table>

<h4>References</h4>

<dl><dt>Idelchik I.E. (1994):</dt>
    <dd><a href=\"http://www.bookfinder.com/dir/i/Handbook_of_Hydraulic_Resistance/0849399084/\"><strong>Handbook
        of Hydraulic Resistance</strong></a>. 3rd edition, Begell House, ISBN
        0-8493-9908-4</dd>
</dl>
</html>"));
        end VesselPortsData;

        connector VesselFluidPorts_b
          "Fluid connector with outlined, large icon to be used for horizontally aligned vectors of FluidPorts (vector dimensions must be added after dragging)"
          extends Interfaces.FluidPort;
          annotation (defaultComponentName="ports_b",
                      Diagram(coordinateSystem(
                preserveAspectRatio=false,
                extent={{-50,-200},{50,200}},
                initialScale=0.2), graphics={
                Text(extent={{-75,130},{75,100}}, textString="%name"),
                Rectangle(
                  extent={{-25,100},{25,-100}}),
                Ellipse(
                  extent={{-22,100},{-10,-100}},
                  fillColor={0,127,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-20,-69},{-12,69}},
                  lineColor={0,127,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-6,100},{6,-100}},
                  fillColor={0,127,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{10,100},{22,-100}},
                  fillColor={0,127,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-4,-69},{4,69}},
                  lineColor={0,127,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{12,-69},{20,69}},
                  lineColor={0,127,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}),
               Icon(coordinateSystem(
                preserveAspectRatio=false,
                extent={{-50,-200},{50,200}},
                initialScale=0.2), graphics={
                Rectangle(
                  extent={{-50,200},{50,-200}},
                  lineColor={0,127,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-44,200},{-20,-200}},
                  fillColor={0,127,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-12,200},{12,-200}},
                  fillColor={0,127,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{20,200},{44,-200}},
                  fillColor={0,127,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-39,-118.5},{-25,113}},
                  lineColor={0,127,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-7,-118.5},{7,113}},
                  lineColor={0,127,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{25,-117.5},{39,114}},
                  lineColor={0,127,255},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}));
        end VesselFluidPorts_b;
      end BaseClasses;
      annotation (Documentation(info="<html>

</html>"));
    end Vessels;

    package Valves "Components for the regulation and control of fluid flow"
        extends Modelica.Icons.VariantsPackage;

      model ValveDiscrete "Valve for water/steam flows with linear pressure drop"
        extends Modelica.Fluid.Interfaces.PartialTwoPortTransport;
        parameter SI.AbsolutePressure dp_nominal
          "Nominal pressure drop at full opening=1"
          annotation(Dialog(group="Nominal operating point"));
        parameter Medium.MassFlowRate m_flow_nominal
          "Nominal mass flowrate at full opening=1";
        final parameter Types.HydraulicConductance k = m_flow_nominal/dp_nominal
          "Hydraulic conductance at full opening=1";
        Modelica.Blocks.Interfaces.BooleanInput open
        annotation (Placement(transformation(
              origin={0,80},
              extent={{-20,-20},{20,20}},
              rotation=270)));
        parameter Real opening_min(min=0)=0
          "Remaining opening if closed, causing small leakage flow";
      equation
        m_flow = if open then 1*k*dp else opening_min*k*dp;

        // Isenthalpic state transformation (no storage and no loss of energy)
        port_a.h_outflow = inStream(port_b.h_outflow);
        port_b.h_outflow = inStream(port_a.h_outflow);

      annotation (
        Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}), graphics={
              Line(points={{0,50},{0,0}}),
              Rectangle(
                extent={{-20,60},{20,50}},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-100,50},{100,-50},{100,50},{0,0},{-100,-50},{-100,50}},
                fillColor=DynamicSelect({255,255,255}, if open then {0,255,0} else {255,255,255}),
                fillPattern=FillPattern.Solid)}),
        Documentation(info="<html>
<p>
This very simple model provides a (small) pressure drop which is proportional to the flowrate if the Boolean open signal is <strong>true</strong>. Otherwise, the mass flow rate is zero. If opening_min > 0, a small leakage mass flow rate occurs when open = <strong>false</strong>.
</p>
<p>This model can be used for simplified modelling of on-off valves, when it is not important to accurately describe the pressure loss when the valve is open. Although the medium model is not used to determine the pressure loss, it must be nevertheless be specified, so that the fluid ports can be connected to other components using the same medium model.</p>
<p>The model is adiabatic (no heat losses to the ambient) and neglects changes in kinetic energy from the inlet to the outlet.</p>
<p>
In a diagram animation, the valve is shown in \"green\", when
it is open.
</p>
</html>", revisions="<html>
<ul>
<li><em>Nov 2005</em>
    by Katja Poschlad (based on ValveLinear).</li>
</ul>
</html>"));
      end ValveDiscrete;
      annotation (Documentation(info="<html>

</html>"));
    end Valves;

    package Sources "Define fixed or prescribed boundary conditions"
      extends Modelica.Icons.SourcesPackage;

      model FixedBoundary "Boundary source component"
        import Modelica.Media.Interfaces.Choices.IndependentVariables;
        extends Sources.BaseClasses.PartialSource;
        parameter Boolean use_p=true "Select p or d"
          annotation (Evaluate = true,
                      Dialog(group = "Boundary pressure or boundary density"));
        parameter Medium.AbsolutePressure p=Medium.p_default "Boundary pressure"
          annotation (Dialog(group = "Boundary pressure or boundary density",
                             enable = use_p));
      parameter Medium.Density d=
       (if use_T then Medium.density_pTX(
                        Medium.p_default,Medium.T_default,Medium.X_default)
        else Medium.density_phX(
                        Medium.p_default,Medium.h_default,Medium.X_default))
          "Boundary density"
          annotation (Dialog(group = "Boundary pressure or boundary density",
                             enable=not use_p));
        parameter Boolean use_T=true "Select T or h"
          annotation (Evaluate = true,
                      Dialog(group = "Boundary temperature or boundary specific enthalpy"));
        parameter Medium.Temperature T = Medium.T_default "Boundary temperature"
          annotation (Dialog(group = "Boundary temperature or boundary specific enthalpy",
                             enable = use_T));
        parameter Medium.SpecificEnthalpy h = Medium.h_default
          "Boundary specific enthalpy"
          annotation (Dialog(group="Boundary temperature or boundary specific enthalpy",
                      enable = not use_T));
        parameter Medium.MassFraction X[Medium.nX](
             quantity=Medium.substanceNames) = Medium.X_default
          "Boundary mass fractions m_i/m"
          annotation (Dialog(group = "Only for multi-substance flow", enable=Medium.nXi > 0));
        parameter Medium.ExtraProperty C[Medium.nC](
             quantity=Medium.extraPropertiesNames) = Medium.C_default
          "Boundary trace substances"
          annotation (Dialog(group = "Only for trace-substance flow", enable=Medium.nC > 0));
    protected
        Medium.ThermodynamicState state;
      equation
        Modelica.Fluid.Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
                                              Medium.singleState, use_p, X,
                                              "FixedBoundary");
        if use_p or Medium.singleState then
           // p given
           if use_T then
              // p,T,X given
              state = Medium.setState_pTX(p, T, X);
           else
              // p,h,X given
              state = Medium.setState_phX(p, h, X);
           end if;

           if Medium.ThermoStates == IndependentVariables.dTX then
              medium.d = Medium.density(state);
           else
              medium.p = Medium.pressure(state);
           end if;

           if Medium.ThermoStates == IndependentVariables.ph or
              Medium.ThermoStates == IndependentVariables.phX then
              medium.h = Medium.specificEnthalpy(state);
           else
              medium.T = Medium.temperature(state);
           end if;

        else
           // d given
           if use_T then
              // d,T,X given
              state = Medium.setState_dTX(d, T, X);

              if Medium.ThermoStates == IndependentVariables.dTX then
                 medium.d = Medium.density(state);
              else
                 medium.p = Medium.pressure(state);
              end if;

              if Medium.ThermoStates == IndependentVariables.ph or
                 Medium.ThermoStates == IndependentVariables.phX then
                 medium.h = Medium.specificEnthalpy(state);
              else
                 medium.T = Medium.temperature(state);
              end if;

           else
              // d,h,X given
              medium.d = d;
              medium.h = h;
              state = Medium.setState_dTX(d,T,X);
           end if;
        end if;

        medium.Xi = X[1:Medium.nXi];

        ports.C_outflow = fill(C, nPorts);
        annotation (defaultComponentName="boundary",
          Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-100,-100},{100,100}}), graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                fillPattern=FillPattern.Sphere,
                fillColor={0,127,255}), Text(
                extent={{-150,110},{150,150}},
                textString="%name",
                textColor={0,0,255})}),
          Documentation(info="<html>
<p>
Model <strong>FixedBoundary</strong> defines constant values for boundary conditions:
</p>
<ul>
<li> Boundary pressure or boundary density.</li>
<li> Boundary temperature or boundary specific enthalpy.</li>
<li> Boundary composition (only for multi-substance or trace-substance flow).</li>
</ul>
<p>
Note, that boundary temperature, density, specific enthalpy,
mass fractions and trace substances have only an effect if the mass flow
is from the Boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary pressure, do not have an effect.
</p>
</html>"));
      end FixedBoundary;

      model MassFlowSource_T
        "Ideal flow source that produces a prescribed mass flow with prescribed temperature, mass fraction and trace substances"
        import Modelica.Media.Interfaces.Choices.IndependentVariables;
        extends Sources.BaseClasses.PartialFlowSource;
        parameter Boolean use_m_flow_in = false
          "Get the mass flow rate from the input connector"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true));
        parameter Boolean use_T_in= false
          "Get the temperature from the input connector"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true));
        parameter Boolean use_X_in = false
          "Get the composition from the input connector"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true));
        parameter Boolean use_C_in = false
          "Get the trace substances from the input connector"
          annotation(Evaluate=true, HideResult=true, choices(checkBox=true));
        parameter Medium.MassFlowRate m_flow = 0
          "Fixed mass flow rate going out of the fluid port"
          annotation (Dialog(enable = not use_m_flow_in));
        parameter Medium.Temperature T = Medium.T_default
          "Fixed value of temperature"
          annotation (Dialog(enable = not use_T_in));
        parameter Medium.MassFraction X[Medium.nX] = Medium.X_default
          "Fixed value of composition"
          annotation (Dialog(enable = (not use_X_in) and Medium.nXi > 0));
        parameter Medium.ExtraProperty C[Medium.nC](
             quantity=Medium.extraPropertiesNames) = Medium.C_default
          "Fixed values of trace substances"
          annotation (Dialog(enable = (not use_C_in) and Medium.nC > 0));
        Modelica.Blocks.Interfaces.RealInput m_flow_in(unit="kg/s") if use_m_flow_in
          "Prescribed mass flow rate"
          annotation (Placement(transformation(extent={{-120,60},{-80,100}}), iconTransformation(extent={{-120,60},{-80,100}})));
        Modelica.Blocks.Interfaces.RealInput T_in(unit="K") if use_T_in
          "Prescribed fluid temperature"
          annotation (Placement(transformation(extent={{-140,20},{-100,60}}), iconTransformation(extent={{-140,20},{-100,60}})));
        Modelica.Blocks.Interfaces.RealInput X_in[Medium.nX](each unit="1") if use_X_in
          "Prescribed fluid composition"
          annotation (Placement(transformation(extent={{-140,-60},{-100,-20}})));
        Modelica.Blocks.Interfaces.RealInput C_in[Medium.nC] if use_C_in
          "Prescribed boundary trace substances"
          annotation (Placement(transformation(extent={{-120,-100},{-80,-60}})));
    protected
        Modelica.Blocks.Interfaces.RealInput m_flow_in_internal(unit="kg/s")
          "Needed to connect to conditional connector";
        Modelica.Blocks.Interfaces.RealInput T_in_internal(unit="K")
          "Needed to connect to conditional connector";
        Modelica.Blocks.Interfaces.RealInput X_in_internal[Medium.nX](each unit="1")
          "Needed to connect to conditional connector";
        Modelica.Blocks.Interfaces.RealInput C_in_internal[Medium.nC]
          "Needed to connect to conditional connector";
      equation
        Utilities.checkBoundary(Medium.mediumName, Medium.substanceNames,
          Medium.singleState, true, X_in_internal, "MassFlowSource_T");
        connect(m_flow_in, m_flow_in_internal);
        connect(T_in, T_in_internal);
        connect(X_in, X_in_internal);
        connect(C_in, C_in_internal);
        if not use_m_flow_in then
          m_flow_in_internal = m_flow;
        end if;
        if not use_T_in then
          T_in_internal = T;
        end if;
        if not use_X_in then
          X_in_internal = X;
        end if;
        if not use_C_in then
          C_in_internal = C;
        end if;
        if Medium.ThermoStates == IndependentVariables.ph or
           Medium.ThermoStates == IndependentVariables.phX then
           medium.h = Medium.specificEnthalpy(Medium.setState_pTX(medium.p, T_in_internal, X_in_internal));
        else
           medium.T = T_in_internal;
        end if;
        sum(ports.m_flow) = -m_flow_in_internal;
        medium.Xi = X_in_internal[1:Medium.nXi];
        ports.C_outflow = fill(C_in_internal, nPorts);
        annotation (defaultComponentName="boundary",
          Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics={
              Rectangle(
                extent={{35,45},{100,-45}},
                fillPattern=FillPattern.HorizontalCylinder,
                fillColor={0,127,255}),
              Ellipse(
                extent={{-100,80},{60,-80}},
                lineColor={0,0,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-60,70},{60,0},{-60,-68},{-60,70}},
                lineColor={0,0,255},
                fillColor={0,0,255},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-54,32},{16,-30}},
                textColor={255,0,0},
                textString="m"),
              Text(
                extent={{-150,130},{150,170}},
                textString="%name",
                textColor={0,0,255}),
              Ellipse(
                extent={{-26,30},{-18,22}},
                lineColor={255,0,0},
                fillColor={255,0,0},
                fillPattern=FillPattern.Solid),
              Text(
                visible=use_m_flow_in,
                extent={{-185,132},{-45,100}},
                textString="m_flow"),
              Text(
                visible=use_T_in,
                extent={{-111,71},{-71,37}},
                textString="T"),
              Text(
                visible=use_X_in,
                extent={{-153,-44},{-33,-72}},
                textString="X"),
              Text(
                visible=use_C_in,
                extent={{-155,-98},{-35,-126}},
                textString="C")}),
          Documentation(info="<html>
<p>
Models an ideal flow source, with prescribed values of flow rate, temperature, composition and trace substances:
</p>
<ul>
<li> Prescribed mass flow rate.</li>
<li> Prescribed temperature.</li>
<li> Boundary composition (only for multi-substance or trace-substance flow).</li>
</ul>
<p>If <code>use_m_flow_in</code> is false (default option), the <code>m_flow</code> parameter
is used as boundary pressure, and the <code>m_flow_in</code> input connector is disabled; if <code>use_m_flow_in</code> is true, then the <code>m_flow</code> parameter is ignored, and the value provided by the input connector is used instead.</p>
<p>The same thing goes for the temperature and composition</p>
<p>
Note, that boundary temperature,
mass fractions and trace substances have only an effect if the mass flow
is from the boundary into the port. If mass is flowing from
the port into the boundary, the boundary definitions,
with exception of boundary flow rate, do not have an effect.
</p>
</html>"));
      end MassFlowSource_T;

      package BaseClasses "Base classes used in the Sources package (only of interest to build new component models)"
        extends Modelica.Icons.BasesPackage;

      partial model PartialSource
          "Partial component source with one fluid connector"
          import Modelica.Constants;

        parameter Integer nPorts=0 "Number of ports" annotation(Dialog(connectorSizing=true));

        replaceable package Medium =
            Modelica.Media.Interfaces.PartialMedium
            "Medium model within the source"
           annotation (choicesAllMatching=true);

        Medium.BaseProperties medium "Medium in the source";

        Interfaces.FluidPorts_b ports[nPorts](
                           redeclare each package Medium = Medium,
                           m_flow(each max=if flowDirection==Types.PortFlowDirection.Leaving then 0 else
                                           +Constants.inf,
                                  each min=if flowDirection==Types.PortFlowDirection.Entering then 0 else
                                           -Constants.inf))
          annotation (Placement(transformation(extent={{90,40},{110,-40}})));
      protected
        parameter Types.PortFlowDirection flowDirection=
                         Types.PortFlowDirection.Bidirectional
            "Allowed flow direction" annotation(Evaluate=true, Dialog(tab="Advanced"));
      equation
        // Only one connection allowed to a port to avoid unwanted ideal mixing
        for i in 1:nPorts loop
          assert(cardinality(ports[i]) <= 1,"
each ports[i] of boundary shall at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts to add an additional port.
");

           ports[i].p          = medium.p;
           ports[i].h_outflow  = medium.h;
           ports[i].Xi_outflow = medium.Xi;
        end for;

        annotation (defaultComponentName="boundary", Documentation(info="<html>
<p>
Partial component to model the <strong>volume interface</strong> of a <strong>source</strong>
component, such as a mass flow source. The essential
features are:
</p>
<ul>
<li> The pressure in the connection port (= ports.p) is identical to the
     pressure in the volume.</li>
<li> The outflow enthalpy rate (= port.h_outflow) and the composition of the
     substances (= port.Xi_outflow) are identical to the respective values in the volume.</li>
</ul>
</html>"));
      end PartialSource;

      partial model PartialFlowSource
          "Partial component source with one fluid connector"
        import Modelica.Constants;

        parameter Integer nPorts=0 "Number of ports" annotation(Dialog(connectorSizing=true));

        replaceable package Medium =
            Modelica.Media.Interfaces.PartialMedium
            "Medium model within the source"
           annotation (choicesAllMatching=true);

        Medium.BaseProperties medium "Medium in the source";

        Interfaces.FluidPort_b ports[nPorts](
                           redeclare each package Medium = Medium,
                           m_flow(each max=if flowDirection==Types.PortFlowDirection.Leaving then 0 else
                                           +Constants.inf,
                                  each min=if flowDirection==Types.PortFlowDirection.Entering then 0 else
                                           -Constants.inf))
          annotation (Placement(transformation(extent={{90,10},{110,-10}})));
      protected
        parameter Types.PortFlowDirection flowDirection=
                         Types.PortFlowDirection.Bidirectional
            "Allowed flow direction" annotation(Evaluate=true, Dialog(tab="Advanced"));
      equation
        assert(abs(sum(abs(ports.m_flow)) - max(abs(ports.m_flow))) <= Modelica.Constants.small, "FlowSource only supports one connection with flow");
        assert(nPorts > 0, "At least one port needs to be present (nPorts > 0), otherwise the model is singular");
        // Only one connection allowed to a port to avoid unwanted ideal mixing
        for i in 1:nPorts loop
          assert(cardinality(ports[i]) <= 1,"
each ports[i] of boundary shall at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts to add an additional port.
");        ports[i].p          = medium.p;
           ports[i].h_outflow  = medium.h;
           ports[i].Xi_outflow = medium.Xi;
        end for;

        annotation (defaultComponentName="boundary", Documentation(info="<html>
<p>
Partial component to model the <strong>volume interface</strong> of a <strong>source</strong>
component, such as a mass flow source. The essential
features are:
</p>
<ul>
<li> The pressure in the connection port (= ports.p) is identical to the
     pressure in the volume.</li>
<li> The outflow enthalpy rate (= port.h_outflow) and the composition of the
     substances (= port.Xi_outflow) are identical to the respective values in the volume.</li>
</ul>
</html>"));
      end PartialFlowSource;
      end BaseClasses;
      annotation (Documentation(info="<html>
<p>
Package <strong>Sources</strong> contains generic sources for fluid connectors
to define fixed or prescribed ambient conditions.
</p>
</html>"));
    end Sources;

    package Sensors "Ideal sensor components to extract signals from a fluid connector"
      extends Modelica.Icons.SensorsPackage;

      model Temperature "Ideal one port temperature sensor"
          extends Sensors.BaseClasses.PartialAbsoluteSensor;

        Modelica.Blocks.Interfaces.RealOutput T(final quantity="ThermodynamicTemperature",
                                                final unit = "K", displayUnit = "degC", min=0)
          "Temperature in port medium"
          annotation (Placement(transformation(extent={{60,-10},{80,10}})));

      equation
        T = Medium.temperature(Medium.setState_phX(port.p, inStream(port.h_outflow), inStream(port.Xi_outflow)));
      annotation (defaultComponentName="temperature",
          Documentation(info="<html>
<p>
This component monitors the temperature of the fluid passing its port.
The sensor is ideal, i.e., it does not influence the fluid.
</p>
</html>"),     Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}}), graphics={
              Line(points={{0,-70},{0,-100}}, color={0,0,127}),
              Ellipse(
                extent={{-20,-98},{20,-60}},
                lineThickness=0.5,
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-12,40},{12,-68}},
                lineColor={191,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-12,40},{-12,80},{-10,86},{-6,88},{0,90},{6,88},{10,86},{
                    12,80},{12,40},{-12,40}},
                lineThickness=0.5),
              Line(
                points={{-12,40},{-12,-64}},
                thickness=0.5),
              Line(
                points={{12,40},{12,-64}},
                thickness=0.5),
              Line(points={{-40,-20},{-12,-20}}),
              Line(points={{-40,20},{-12,20}}),
              Line(points={{-40,60},{-12,60}}),
              Line(points={{12,0},{60,0}}, color={0,0,127})}),
          Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Ellipse(
                extent={{-20,-88},{20,-50}},
                lineThickness=0.5,
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-12,50},{12,-58}},
                lineColor={191,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-12,50},{-12,90},{-10,96},{-6,98},{0,100},{6,98},{10,96},{
                    12,90},{12,50},{-12,50}},
                lineThickness=0.5),
              Line(
                points={{-12,50},{-12,-54}},
                thickness=0.5),
              Line(
                points={{12,50},{12,-54}},
                thickness=0.5),
              Line(points={{-40,-10},{-12,-10}}),
              Line(points={{-40,30},{-12,30}}),
              Line(points={{-40,70},{-12,70}}),
              Text(
                extent={{126,-30},{6,-60}},
                textString="T"),
              Text(
                extent={{-150,110},{150,150}},
                textString="%name",
                textColor={0,0,255}),
              Line(points={{12,0},{60,0}}, color={0,0,127})}));
      end Temperature;

      model TemperatureTwoPort "Ideal two port temperature sensor"
        extends Sensors.BaseClasses.PartialFlowSensor;

        Modelica.Blocks.Interfaces.RealOutput T( final quantity="ThermodynamicTemperature",
                                                 final unit="K",
                                                 min = 0,
                                                 displayUnit="degC")
          "Temperature of the passing fluid"
          annotation (Placement(transformation(
              origin={0,110},
              extent={{10,-10},{-10,10}},
              rotation=270)));

    protected
        Medium.Temperature T_a_inflow "Temperature of inflowing fluid at port_a";
        Medium.Temperature T_b_inflow
          "Temperature of inflowing fluid at port_b or T_a_inflow, if uni-directional flow";
      equation
        if allowFlowReversal then
           T_a_inflow = Medium.temperature(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
           T_b_inflow = Medium.temperature(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow));
           T = Modelica.Fluid.Utilities.regStep(port_a.m_flow, T_a_inflow, T_b_inflow, m_flow_small);
        else
           T = Medium.temperature(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow));
           T_a_inflow = T;
           T_b_inflow = T;
        end if;
      annotation (defaultComponentName="temperature",
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Line(points={{0,100},{0,50}}, color={0,0,127}),
              Line(points={{-92,0},{100,0}}, color={0,128,255}),
              Ellipse(
                extent={{-20,-68},{20,-30}},
                lineThickness=0.5,
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-12,50},{12,-34}},
                lineColor={191,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Polygon(
                points={{-12,50},{-12,70},{-10,76},{-6,78},{0,80},{6,78},{10,76},{
                    12,70},{12,50},{-12,50}},
                lineThickness=0.5),
              Line(
                points={{-12,50},{-12,-35}},
                thickness=0.5),
              Line(
                points={{12,50},{12,-34}},
                thickness=0.5),
              Line(points={{-40,-10},{-12,-10}}),
              Line(points={{-40,20},{-12,20}}),
              Line(points={{-40,50},{-12,50}}),
              Text(
                extent={{94,122},{0,92}},
                textString="T")}),
        Documentation(info="<html>
<p>
This component monitors the temperature of the passing fluid.
The sensor is ideal, i.e., it does not influence the fluid.
</p>
</html>"));
      end TemperatureTwoPort;

      model SpecificEnthalpyTwoPort
        "Ideal two port sensor for the specific enthalpy"
        extends Sensors.BaseClasses.PartialFlowSensor;
        extends Modelica.Icons.RoundSensor;
        Modelica.Blocks.Interfaces.RealOutput h_out(final quantity="SpecificEnergy",
                                                    final unit="J/kg")
          "Specific enthalpy of the passing fluid"
          annotation (Placement(transformation(
              origin={0,110},
              extent={{10,-10},{-10,10}},
              rotation=270)));

      equation
        if allowFlowReversal then
           h_out = Modelica.Fluid.Utilities.regStep(port_a.m_flow, port_b.h_outflow, port_a.h_outflow, m_flow_small);
        else
           h_out = port_b.h_outflow;
        end if;
      annotation (defaultComponentName="specificEnthalpy",
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Text(
                extent={{102,120},{0,90}},
                textString="h"),
              Line(points={{0,100},{0,70}}, color={0,0,127}),
              Line(points={{-100,0},{-70,0}}, color={0,128,255}),
              Line(points={{70,0},{100,0}}, color={0,128,255})}),
        Documentation(info="<html>
<p>
This component monitors the specific enthalpy of a passing fluid.
The sensor is ideal, i.e., it does not influence the fluid.
</p>
</html>"));
      end SpecificEnthalpyTwoPort;

      model MassFlowRate "Ideal sensor for mass flow rate"
        extends Sensors.BaseClasses.PartialFlowSensor;
        extends Modelica.Icons.RoundSensor;
        Modelica.Blocks.Interfaces.RealOutput m_flow(quantity="MassFlowRate",
                                                     final unit="kg/s")
          "Mass flow rate from port_a to port_b" annotation (Placement(
              transformation(
              origin={0,110},
              extent={{10,-10},{-10,10}},
              rotation=270)));

      equation
        m_flow = port_a.m_flow;
      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
                  100}}), graphics={
              Line(points={{70,0},{100,0}}, color={0,128,255}),
              Text(
                extent={{162,120},{2,90}},
                textString="m_flow"),
              Line(points={{0,100},{0,70}}, color={0,0,127}),
              Line(points={{-100,0},{-70,0}}, color={0,128,255})}),
        Documentation(info="<html>
<p>
This component monitors the mass flow rate flowing from port_a to port_b.
The sensor is ideal, i.e., it does not influence the fluid.
</p>
</html>"));
      end MassFlowRate;

      package BaseClasses "Base classes used in the Sensors package (only of interest to build new component models)"
        extends Modelica.Icons.BasesPackage;

        partial model PartialAbsoluteSensor
          "Partial component to model a sensor that measures a potential variable"

          replaceable package Medium=Modelica.Media.Interfaces.PartialMedium
            "Medium in the sensor"
            annotation(choicesAllMatching=true);

          Modelica.Fluid.Interfaces.FluidPort_a port(redeclare package Medium=Medium, m_flow(min=0))
            annotation (Placement(transformation(
                origin={0,-100},
                extent={{-10,-10},{10,10}},
                rotation=90)));

        equation
          port.m_flow = 0;
          port.h_outflow = Medium.h_default;
          port.Xi_outflow = Medium.X_default[1:Medium.nXi];
          port.C_outflow = zeros(Medium.nC);
          annotation (Documentation(info="<html>
<p>
Partial component to model an <strong>absolute sensor</strong>. Can be used for pressure sensor models.
Use for other properties such as temperature or density is discouraged, because the enthalpy at the connector can have different meanings, depending on the connection topology. Use <code>PartialFlowSensor</code> instead.
as signal.
</p>
</html>"));
        end PartialAbsoluteSensor;

        partial model PartialFlowSensor
          "Partial component to model sensors that measure flow properties"
          extends Modelica.Fluid.Interfaces.PartialTwoPort;

          parameter Medium.MassFlowRate m_flow_nominal = system.m_flow_nominal
            "Nominal value of m_flow = port_a.m_flow"
            annotation(Dialog(tab = "Advanced"));
          parameter Medium.MassFlowRate m_flow_small(min=0) = if system.use_eps_Re then system.eps_m_flow*m_flow_nominal else system.m_flow_small
            "Regularization for bi-directional flow in the region |m_flow| < m_flow_small (m_flow_small > 0 required)"
            annotation(Dialog(tab="Advanced"));

        equation
          // mass balance
          0 = port_a.m_flow + port_b.m_flow;

          // momentum equation (no pressure loss)
          port_a.p = port_b.p;

          // isenthalpic state transformation (no storage and no loss of energy)
          port_a.h_outflow = inStream(port_b.h_outflow);
          port_b.h_outflow = inStream(port_a.h_outflow);

          port_a.Xi_outflow = inStream(port_b.Xi_outflow);
          port_b.Xi_outflow = inStream(port_a.Xi_outflow);

          port_a.C_outflow = inStream(port_b.C_outflow);
          port_b.C_outflow = inStream(port_a.C_outflow);
          annotation (Documentation(info="<html>
<p>
Partial component to model a <strong>sensor</strong> that measures any intensive properties
of a flow, e.g., to get temperature or density in the flow
between fluid connectors.<br>
The model includes zero-volume balance equations. Sensor models inheriting from
this partial class should add a medium instance to calculate the measured property.
</p>
</html>"));
        end PartialFlowSensor;
      end BaseClasses;
      annotation (preferredView="info", Documentation(info="<html>
<p>
Package <strong>Sensors</strong> consists of idealized sensor components that
provide variables of a medium model and/or fluid ports as
output signals. These signals can be, e.g., further processed
with components of the Modelica.Blocks library.
Also more realistic sensor models can be built, by further
processing (e.g., by attaching block Modelica.Blocks.FirstOrder to
model the time constant of the sensor).
</p>

<p>For the thermodynamic state variables temperature, specific enthalpy, specific entropy and density
the fluid library provides two different types of sensors: <strong>regular one port</strong> and <strong>two port</strong> sensors.</p>

<ul>
<li>The <strong>regular one port</strong> sensors have the advantage of easy introduction and removal from a model, as no connections have to be broken.
A potential drawback is that the obtained value jumps as flow reverts.
</li>

<li>The <strong>two port</strong> sensors offer the advantages of an adjustable regularized step function around zero flow.
Moreover the obtained result is restricted to the value flowing into port_a if allowFlowReversal is false.</li>
</ul>

<p>
<a href=\"modelica://Modelica.Fluid.Examples.Explanatory.MeasuringTemperature\">Modelica.Fluid.Examples.Explanatory.MeasuringTemperature</a>
demonstrates the differences between one- and two-port sensor at hand of a
simple example.
</p>
</html>",     revisions="<html>
<ul>
<li><em>22 Dec 2008</em>
    by R;uumldiger Franke<br>
    <ul>
    <li>flow sensors based on Interfaces.PartialTwoPort</li>
    <li>adapted docu to stream connectors, i.e., less need for two port sensors</li>
    </ul>
    </li>
<li><em>4 Dec 2008</em>
    by Michael Wetter<br>
       included sensors for trace substance</li>
<li><em>31 Oct 2007</em>
    by Carsten Heinrich<br>
       updated sensor models, included one and two port sensors for thermodynamic state variables</li>
</ul>
</html>"));
    end Sensors;

    package Interfaces "Interfaces for steady state and unsteady, mixed-phase, multi-substance, incompressible and compressible flow"
      extends Modelica.Icons.InterfacesPackage;

      connector FluidPort
        "Interface for quasi one-dimensional fluid flow in a piping network (incompressible or compressible, one or more phases, one or more substances)"

        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium
          "Medium model" annotation (choicesAllMatching=true);

        flow Medium.MassFlowRate m_flow
          "Mass flow rate from the connection point into the component";
        Medium.AbsolutePressure p "Thermodynamic pressure in the connection point";
        stream Medium.SpecificEnthalpy h_outflow
          "Specific thermodynamic enthalpy close to the connection point if m_flow < 0";
        stream Medium.MassFraction Xi_outflow[Medium.nXi]
          "Independent mixture mass fractions m_i/m close to the connection point if m_flow < 0";
        stream Medium.ExtraProperty C_outflow[Medium.nC]
          "Properties c_i/m close to the connection point if m_flow < 0";
      end FluidPort;

      connector FluidPort_a "Generic fluid connector at design inlet"
        extends FluidPort;
        annotation (defaultComponentName="port_a",
                    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics={Ellipse(
                extent={{-40,40},{40,-40}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid), Text(extent={{-150,110},{150,50}},
                  textString="%name")}),
             Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}}), graphics={Ellipse(
                extent={{-100,100},{100,-100}},
                lineColor={0,127,255},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid), Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid)}));
      end FluidPort_a;

      connector FluidPort_b "Generic fluid connector at design outlet"
        extends FluidPort;
        annotation (defaultComponentName="port_b",
                    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                  -100},{100,100}}), graphics={
              Ellipse(
                extent={{-40,40},{40,-40}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-30,30},{30,-30}},
                lineColor={0,127,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Text(extent={{-150,110},{150,50}}, textString="%name")}),
             Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
                  100,100}}), graphics={
              Ellipse(
                extent={{-100,100},{100,-100}},
                lineColor={0,127,255},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-100,100},{100,-100}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-80,80},{80,-80}},
                lineColor={0,127,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}));
      end FluidPort_b;

      connector FluidPorts_b
        "Fluid connector with outlined, large icon to be used for vectors of FluidPorts (vector dimensions must be added after dragging)"
        extends FluidPort;
        annotation (defaultComponentName="ports_b",
                    Diagram(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-50,-200},{50,200}},
              initialScale=0.2), graphics={
              Text(extent={{-75,130},{75,100}}, textString="%name"),
              Rectangle(
                extent={{-25,100},{25,-100}}),
              Ellipse(
                extent={{-25,90},{25,40}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-25,25},{25,-25}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-25,-40},{25,-90}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-15,-50},{15,-80}},
                lineColor={0,127,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-15,15},{15,-15}},
                lineColor={0,127,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-15,50},{15,80}},
                lineColor={0,127,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}),
             Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-50,-200},{50,200}},
              initialScale=0.2), graphics={
              Rectangle(
                extent={{-50,200},{50,-200}},
                lineColor={0,127,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-50,180},{50,80}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-50,50},{50,-50}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-50,-80},{50,-180}},
                fillColor={0,127,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-30,30},{30,-30}},
                lineColor={0,127,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-30,100},{30,160}},
                lineColor={0,127,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Ellipse(
                extent={{-30,-100},{30,-160}},
                lineColor={0,127,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}));
      end FluidPorts_b;

      partial model PartialTwoPort "Partial component with two ports"
        import Modelica.Constants;
        outer Modelica.Fluid.System system "System wide properties";

        replaceable package Medium =
            Modelica.Media.Interfaces.PartialMedium "Medium in the component"
            annotation (choicesAllMatching = true);

        parameter Boolean allowFlowReversal = system.allowFlowReversal
          "= true to allow flow reversal, false restricts to design direction (port_a -> port_b)"
          annotation(Dialog(tab="Assumptions"), Evaluate=true);

        Modelica.Fluid.Interfaces.FluidPort_a port_a(
                                      redeclare package Medium = Medium,
                           m_flow(min=if allowFlowReversal then -Constants.inf else 0))
          "Fluid connector a (positive design flow direction is from port_a to port_b)"
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
        Modelica.Fluid.Interfaces.FluidPort_b port_b(
                                      redeclare package Medium = Medium,
                           m_flow(max=if allowFlowReversal then +Constants.inf else 0))
          "Fluid connector b (positive design flow direction is from port_a to port_b)"
          annotation (Placement(transformation(extent={{110,-10},{90,10}}), iconTransformation(extent={{110,-10},{90,10}})));
        // Model structure, e.g., used for visualization
    protected
        parameter Boolean port_a_exposesState = false
          "= true if port_a exposes the state of a fluid volume";
        parameter Boolean port_b_exposesState = false
          "= true if port_b.p exposes the state of a fluid volume";
        parameter Boolean showDesignFlowDirection = true
          "= false to hide the arrow in the model icon";

        annotation (
          Documentation(info="<html>
<p>
This partial model defines an interface for components with two ports.
The treatment of the design flow direction and of flow reversal are predefined based on the parameter <code><strong>allowFlowReversal</strong></code>.
The component may transport fluid and may have internal storage for a given fluid <code><strong>Medium</strong></code>.
</p>
<p>
An extending model providing direct access to internal storage of mass or energy through port_a or port_b
should redefine the protected parameters <code><strong>port_a_exposesState</strong></code> and <code><strong>port_b_exposesState</strong></code> appropriately.
This will be visualized at the port icons, in order to improve the understanding of fluid model diagrams.
</p>
</html>"),Icon(coordinateSystem(
                preserveAspectRatio=true,
                extent={{-100,-100},{100,100}}), graphics={
              Polygon(
                points={{20,-70},{60,-85},{20,-100},{20,-70}},
                lineColor={0,128,255},
                fillColor={0,128,255},
                fillPattern=FillPattern.Solid,
                visible=showDesignFlowDirection),
              Polygon(
                points={{20,-75},{50,-85},{20,-95},{20,-75}},
                lineColor={255,255,255},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                visible=allowFlowReversal),
              Line(
                points={{55,-85},{-60,-85}},
                color={0,128,255},
                visible=showDesignFlowDirection),
              Text(
                extent={{-149,-114},{151,-154}},
                textColor={0,0,255},
                textString="%name"),
              Ellipse(
                extent={{-110,26},{-90,-24}},
                fillPattern=FillPattern.Solid,
                visible=port_a_exposesState),
              Ellipse(
                extent={{90,25},{110,-25}},
                fillPattern=FillPattern.Solid,
                visible=port_b_exposesState)}));
      end PartialTwoPort;

    partial model PartialTwoPortTransport
        "Partial element transporting fluid between two ports without storage of mass or energy"

      extends PartialTwoPort(
        final port_a_exposesState=false,
        final port_b_exposesState=false);

      // Advanced
      // Note: value of dp_start shall be refined by derived model, basing on local dp_nominal
      parameter Medium.AbsolutePressure dp_start(min=-Modelica.Constants.inf) = 0.01*system.p_start
          "Guess value of dp = port_a.p - port_b.p"
        annotation(Dialog(tab = "Advanced"));
      parameter Medium.MassFlowRate m_flow_start = system.m_flow_start
          "Guess value of m_flow = port_a.m_flow"
        annotation(Dialog(tab = "Advanced"));
      // Note: value of m_flow_small shall be refined by derived model, basing on local m_flow_nominal
      parameter Medium.MassFlowRate m_flow_small = if system.use_eps_Re then system.eps_m_flow*system.m_flow_nominal else system.m_flow_small
          "Small mass flow rate for regularization of zero flow"
        annotation(Dialog(tab = "Advanced"));

      // Diagnostics
      parameter Boolean show_T = true
          "= true, if temperatures at port_a and port_b are computed"
        annotation(Dialog(tab="Advanced",group="Diagnostics"));
      parameter Boolean show_V_flow = true
          "= true, if volume flow rate at inflowing port is computed"
        annotation(Dialog(tab="Advanced",group="Diagnostics"));

      // Variables
      Medium.MassFlowRate m_flow(
         min=if allowFlowReversal then -Modelica.Constants.inf else 0,
         start = m_flow_start) "Mass flow rate in design flow direction";
      SI.Pressure dp(start=dp_start)
          "Pressure difference between port_a and port_b (= port_a.p - port_b.p)";

      SI.VolumeFlowRate V_flow=
          m_flow/Modelica.Fluid.Utilities.regStep(m_flow,
                      Medium.density(state_a),
                      Medium.density(state_b),
                      m_flow_small) if show_V_flow
          "Volume flow rate at inflowing port (positive when flow from port_a to port_b)";

      Medium.Temperature port_a_T=
          Modelica.Fluid.Utilities.regStep(port_a.m_flow,
                      Medium.temperature(state_a),
                      Medium.temperature(Medium.setState_phX(port_a.p, port_a.h_outflow, port_a.Xi_outflow)),
                      m_flow_small) if show_T
          "Temperature close to port_a, if show_T = true";
      Medium.Temperature port_b_T=
          Modelica.Fluid.Utilities.regStep(port_b.m_flow,
                      Medium.temperature(state_b),
                      Medium.temperature(Medium.setState_phX(port_b.p, port_b.h_outflow, port_b.Xi_outflow)),
                      m_flow_small) if show_T
          "Temperature close to port_b, if show_T = true";
    protected
      Medium.ThermodynamicState state_a "State for medium inflowing through port_a";
      Medium.ThermodynamicState state_b "State for medium inflowing through port_b";
    equation
      // medium states
      state_a = Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow));
      state_b = Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow));

      // Pressure drop in design flow direction
      dp = port_a.p - port_b.p;

      // Design direction of mass flow rate
      m_flow = port_a.m_flow;
      assert(m_flow > -m_flow_small or allowFlowReversal, "Reversing flow occurs even though allowFlowReversal is false");

      // Mass balance (no storage)
      port_a.m_flow + port_b.m_flow = 0;

      // Transport of substances
      port_a.Xi_outflow = inStream(port_b.Xi_outflow);
      port_b.Xi_outflow = inStream(port_a.Xi_outflow);

      port_a.C_outflow = inStream(port_b.C_outflow);
      port_b.C_outflow = inStream(port_a.C_outflow);

      annotation (
        Documentation(info="<html>
<p>
This component transports fluid between its two ports, without storing mass or energy.
Energy may be exchanged with the environment though, e.g., in the form of work.
<code>PartialTwoPortTransport</code> is intended as base class for devices like orifices, valves and simple fluid machines.</p>
<p>
Three equations need to be added by an extending class using this component:
</p>
<ul>
<li>the momentum balance specifying the relationship between the pressure drop <code>dp</code> and the mass flow rate <code>m_flow</code>,</li>
<li><code>port_b.h_outflow</code> for flow in design direction, and</li>
<li><code>port_a.h_outflow</code> for flow in reverse direction.</li>
</ul>
<p>
Moreover appropriate values shall be assigned to the following parameters:
</p>
<ul>
<li><code>dp_start</code> for a guess of the pressure drop</li>
<li><code>m_flow_small</code> for regularization of zero flow.</li>
</ul>
</html>"));
    end PartialTwoPortTransport;

      connector HeatPorts_a
        "HeatPort connector with filled, large icon to be used for vectors of HeatPorts (vector dimensions must be added after dragging)"
        extends Modelica.Thermal.HeatTransfer.Interfaces.HeatPort;
        annotation (defaultComponentName="heatPorts_a",
             Icon(coordinateSystem(
              preserveAspectRatio=false,
              extent={{-200,-50},{200,50}},
              initialScale=0.2), graphics={
              Rectangle(
                extent={{-201,50},{200,-50}},
                lineColor={127,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-171,45},{-83,-45}},
                lineColor={127,0,0},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{-45,45},{43,-45}},
                lineColor={127,0,0},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid),
              Rectangle(
                extent={{82,45},{170,-45}},
                lineColor={127,0,0},
                fillColor={127,0,0},
                fillPattern=FillPattern.Solid)}));
      end HeatPorts_a;

      partial model PartialHeatTransfer "Common interface for heat transfer models"

        // Parameters
        replaceable package Medium=Modelica.Media.Interfaces.PartialMedium
          "Medium in the component"
          annotation(Dialog(tab="Internal Interface",enable=false));

        parameter Integer n=1 "Number of heat transfer segments"
          annotation(Dialog(tab="Internal Interface",enable=false), Evaluate=true);

        // Inputs provided to heat transfer model
        input Medium.ThermodynamicState[n] states
          "Thermodynamic states of flow segments";

        input SI.Area[n] surfaceAreas "Heat transfer areas";

        // Outputs defined by heat transfer model
        output SI.HeatFlowRate[n] Q_flows "Heat flow rates";

        // Parameters
        parameter Boolean use_k = false
          "= true to use k value for thermal isolation"
          annotation(Dialog(tab="Internal Interface",enable=false));
        parameter SI.CoefficientOfHeatTransfer k = 0
          "Heat transfer coefficient to ambient"
          annotation(Dialog(group="Ambient"),Evaluate=true);
        parameter SI.Temperature T_ambient = system.T_ambient "Ambient temperature"
          annotation(Dialog(group="Ambient"));
        outer Modelica.Fluid.System system "System wide properties";

        // Heat ports
        Modelica.Fluid.Interfaces.HeatPorts_a[n] heatPorts
          "Heat port to component boundary"
          annotation (Placement(transformation(extent={{-10,60},{10,80}}), iconTransformation(extent={{-20,60},{20,80}})));

        // Variables
        SI.Temperature[n] Ts = Medium.temperature(states)
          "Temperatures defined by fluid states";

      equation
        if use_k then
          Q_flows = heatPorts.Q_flow + {k*surfaceAreas[i]*(T_ambient - heatPorts[i].T) for i in 1:n};
        else
          Q_flows = heatPorts.Q_flow;
        end if;

        annotation (Documentation(info="<html>
<p>
This component is a common interface for heat transfer models. The heat flow rates <code>Q_flows[n]</code> through the boundaries of n flow segments
are obtained as function of the thermodynamic <code>states</code> of the flow segments for a given fluid <code>Medium</code>,
the <code>surfaceAreas[n]</code> and the boundary temperatures <code>heatPorts[n].T</code>.
</p>
<p>
The heat loss coefficient <code>k</code> can be used to model a thermal isolation between <code>heatPorts.T</code> and <code>T_ambient</code>.</p>
<p>
An extending model implementing this interface needs to define one equation: the relation between the predefined fluid temperatures <code>Ts[n]</code>,
the boundary temperatures <code>heatPorts[n].T</code>, and the heat flow rates <code>Q_flows[n]</code>.
</p>
</html>"));
      end PartialHeatTransfer;

        partial model PartialLumpedVolume
        "Lumped volume with mass and energy balance"
        import Modelica.Fluid.Types;
        import Modelica.Fluid.Types.Dynamics;
        import Modelica.Media.Interfaces.Choices.IndependentVariables;

          outer Modelica.Fluid.System system "System properties";
          replaceable package Medium =
            Modelica.Media.Interfaces.PartialMedium "Medium in the component"
              annotation (choicesAllMatching = true);

          // Inputs provided to the volume model
          input SI.Volume fluidVolume "Volume";

          // Assumptions
          parameter Types.Dynamics energyDynamics=system.energyDynamics
          "Formulation of energy balance"
            annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
          parameter Types.Dynamics massDynamics=system.massDynamics
          "Formulation of mass balance"
            annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
          final parameter Types.Dynamics substanceDynamics=massDynamics
          "Formulation of substance balance"
            annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
          final parameter Types.Dynamics traceDynamics=massDynamics
          "Formulation of trace substance balance"
            annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

          // Initialization
          parameter Medium.AbsolutePressure p_start = system.p_start
          "Start value of pressure"
            annotation(Dialog(tab = "Initialization"));
          parameter Boolean use_T_start = true
          "= true, use T_start, otherwise h_start"
            annotation(Dialog(tab = "Initialization"), Evaluate=true);
          parameter Medium.Temperature T_start=
            if use_T_start then system.T_start else Medium.temperature_phX(p_start,h_start,X_start)
          "Start value of temperature"
            annotation(Dialog(tab = "Initialization", enable = use_T_start));
          parameter Medium.SpecificEnthalpy h_start=
            if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default
          "Start value of specific enthalpy"
            annotation(Dialog(tab = "Initialization", enable = not use_T_start));
          parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default
          "Start value of mass fractions m_i/m"
            annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
          parameter Medium.ExtraProperty C_start[Medium.nC](
               quantity=Medium.extraPropertiesNames) = Medium.C_default
          "Start value of trace substances"
            annotation (Dialog(tab="Initialization", enable=Medium.nC > 0));

          Medium.BaseProperties medium(
            preferredMediumStates = (if energyDynamics == Dynamics.SteadyState and
                                        massDynamics   == Dynamics.SteadyState then false else true),
            p(start=p_start),
            h(start=h_start),
            T(start=T_start),
            Xi(start=X_start[1:Medium.nXi]));
          SI.Energy U "Internal energy of fluid";
          SI.Mass m "Mass of fluid";
          SI.Mass[Medium.nXi] mXi "Masses of independent components in the fluid";
          SI.Mass[Medium.nC] mC "Masses of trace substances in the fluid";
          // C need to be added here because unlike for Xi, which has medium.Xi,
          // there is no variable medium.C
          Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";

          // variables that need to be defined by an extending class
          SI.MassFlowRate mb_flow "Mass flows across boundaries";
          SI.MassFlowRate[Medium.nXi] mbXi_flow
          "Substance mass flows across boundaries";
          Medium.ExtraPropertyFlowRate[Medium.nC] mbC_flow
          "Trace substance mass flows across boundaries";
          SI.EnthalpyFlowRate Hb_flow
          "Enthalpy flow across boundaries or energy source/sink";
          SI.HeatFlowRate Qb_flow
          "Heat flow across boundaries or energy source/sink";
          SI.Power Wb_flow "Work flow across boundaries or source term";
    protected
          parameter Boolean initialize_p = not Medium.singleState
          "= true to set up initial equations for pressure";
          Real[Medium.nC] mC_scaled(min=fill(Modelica.Constants.eps, Medium.nC))
          "Scaled masses of trace substances in the fluid";
        equation
          assert(not (energyDynamics<>Dynamics.SteadyState and massDynamics==Dynamics.SteadyState) or Medium.singleState,
                 "Bad combination of dynamics options and Medium not conserving mass if fluidVolume is fixed.");

          // Total quantities
          m = fluidVolume*medium.d;
          mXi = m*medium.Xi;
          U = m*medium.u;
          mC = m*C;

          // Energy and mass balances
          if energyDynamics == Dynamics.SteadyState then
            0 = Hb_flow + Qb_flow + Wb_flow;
          else
            der(U) = Hb_flow + Qb_flow + Wb_flow;
          end if;

          if massDynamics == Dynamics.SteadyState then
            0 = mb_flow;
          else
            der(m) = mb_flow;
          end if;

          if substanceDynamics == Dynamics.SteadyState then
            zeros(Medium.nXi) = mbXi_flow;
          else
            der(mXi) = mbXi_flow;
          end if;

          if traceDynamics == Dynamics.SteadyState then
            zeros(Medium.nC)  = mbC_flow;
          else
            der(mC_scaled) = mbC_flow./Medium.C_nominal;
          end if;
            mC = mC_scaled.*Medium.C_nominal;

        initial equation
          // initialization of balances
          if energyDynamics == Dynamics.FixedInitial then
            /*
    if use_T_start then
      medium.T = T_start;
    else
      medium.h = h_start;
    end if;
    */
            if Medium.ThermoStates == IndependentVariables.ph or
               Medium.ThermoStates == IndependentVariables.phX then
               medium.h = h_start;
            else
               medium.T = T_start;
            end if;
          elseif energyDynamics == Dynamics.SteadyStateInitial then
            /*
    if use_T_start then
      der(medium.T) = 0;
    else
      der(medium.h) = 0;
    end if;
    */
            if Medium.ThermoStates == IndependentVariables.ph or
               Medium.ThermoStates == IndependentVariables.phX then
               der(medium.h) = 0;
            else
               der(medium.T) = 0;
            end if;
          end if;

          if massDynamics == Dynamics.FixedInitial then
            if initialize_p then
              medium.p = p_start;
            end if;
          elseif massDynamics == Dynamics.SteadyStateInitial then
            if initialize_p then
              der(medium.p) = 0;
            end if;
          end if;

          if substanceDynamics == Dynamics.FixedInitial then
            medium.Xi = X_start[1:Medium.nXi];
          elseif substanceDynamics == Dynamics.SteadyStateInitial then
            der(medium.Xi) = zeros(Medium.nXi);
          end if;

          if traceDynamics == Dynamics.FixedInitial then
            mC_scaled = m*C_start[1:Medium.nC]./Medium.C_nominal;
          elseif traceDynamics == Dynamics.SteadyStateInitial then
            der(mC_scaled) = zeros(Medium.nC);
          end if;

          annotation (
            Documentation(info="<html>
<p>
Interface and base class for an ideally mixed fluid volume with the ability to store mass and energy.
The following boundary flow and source terms are part of the energy balance and must be specified in an extending class:
</p>
<ul>
<li><code><strong>Qb_flow</strong></code>, e.g., convective or latent heat flow rate across segment boundary, and</li>
<li><code><strong>Wb_flow</strong></code>, work term, e.g., p*der(fluidVolume) if the volume is not constant.</li>
</ul>
<p>
The component volume <code><strong>fluidVolume</strong></code> is an input that needs to be set in the extending class to complete the model.
</p>
<p>
Further source terms must be defined by an extending class for fluid flow across the segment boundary:
</p>
<ul>
<li><code><strong>Hb_flow</strong></code>, enthalpy flow,</li>
<li><code><strong>mb_flow</strong></code>, mass flow,</li>
<li><code><strong>mbXi_flow</strong></code>, substance mass flow, and</li>
<li><code><strong>mbC_flow</strong></code>, trace substance mass flow.</li>
</ul>
</html>"));
        end PartialLumpedVolume;
      annotation (Documentation(info="<html>

</html>",     revisions="<html>
<ul>
<li><em>June 9th, 2008</em>
       by Michael Sielemann: Introduced stream keyword after decision at 57th Design Meeting (Lund).</li>
<li><em>May 30, 2007</em>
       by Christoph Richter: moved everything back to its original position in Modelica.Fluid.</li>
<li><em>Apr. 20, 2007</em>
       by Christoph Richter: moved parts of the original package from Modelica.Fluid
       to the development branch of Modelica 2.2.2.</li>
<li><em>Nov. 2, 2005</em>
       by Francesco Casella: restructured after 45th Design Meeting.</li>
<li><em>Nov. 20-21, 2002</em>
       by Hilding Elmqvist, Mike Tiller, Allan Watson, John Batteh, Chuck Newman,
       Jonas Eborn: Improved at the 32nd Modelica Design Meeting.</li>
<li><em>Nov. 11, 2002</em>
       by Hilding Elmqvist, Martin Otter: improved version.</li>
<li><em>Nov. 6, 2002</em>
       by Hilding Elmqvist: first version.</li>
<li><em>Aug. 11, 2002</em>
       by Martin Otter: Improved according to discussion with Hilding
       Elmqvist and Hubertus Tummescheit.<br>
       The PortVicinity model is manually
       expanded in the base models.<br>
       The Volume used for components is renamed
       PartialComponentVolume.<br>
       A new volume model \"Fluid.Components.PortVolume\"
       introduced that has the medium properties of the port to which it is
       connected.<br>
       Fluid.Interfaces.PartialTwoPortTransport is a component
       for elementary two port transport elements, whereas PartialTwoPort
       is a component for a container component.</li>
</ul>
</html>"));
    end Interfaces;

    package Types "Common types for fluid models"
      extends Modelica.Icons.TypesPackage;

      type HydraulicConductance = Modelica.Icons.TypeReal (
          final quantity="HydraulicConductance",
          final unit="kg/(s.Pa)") "Real type for hydraulic conductance";

      type Dynamics = enumeration(
          DynamicFreeInitial
            "DynamicFreeInitial -- Dynamic balance, Initial guess value",
          FixedInitial "FixedInitial -- Dynamic balance, Initial value fixed",
          SteadyStateInitial
            "SteadyStateInitial -- Dynamic balance, Steady state initial with guess value",
          SteadyState "SteadyState -- Steady state balance, Initial guess value")
        "Enumeration to define definition of balance equations"
      annotation (Documentation(info="<html>
<p>
Enumeration to define the formulation of balance equations
(to be selected via choices menu):
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th><strong>Dynamics.</strong></th><th><strong>Meaning</strong></th></tr>
<tr><td>DynamicFreeInitial</td><td>Dynamic balance, Initial guess value</td></tr>

<tr><td>FixedInitial</td><td>Dynamic balance, Initial value fixed</td></tr>

<tr><td>SteadyStateInitial</td><td>Dynamic balance, Steady state initial with guess value</td></tr>

<tr><td>SteadyState</td><td>Steady state balance, Initial guess value</td></tr>
</table>

<p>
The enumeration \"Dynamics\" is used for the mass, energy and momentum balance equations
respectively. The exact meaning for the three balance equations is stated in the following
tables:
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><td colspan=\"3\"><strong>Mass balance</strong> </td></tr>
<tr><td><strong>Dynamics.</strong></td>
    <td><strong>Balance equation</strong></td>
    <td><strong>Initial condition</strong></td></tr>

<tr><td> DynamicFreeInitial</td>
    <td> no restrictions </td>
    <td> no initial conditions </td></tr>

<tr><td> FixedInitial</td>
    <td> no restrictions </td>
    <td> <strong>if</strong> Medium.singleState <strong>then</strong><br>
         &nbsp;&nbsp;no initial condition<br>
         <strong>else</strong> p=p_start </td></tr>

<tr><td> SteadyStateInitial</td>
    <td> no restrictions </td>
    <td> <strong>if</strong> Medium.singleState <strong>then</strong><br>
         &nbsp;&nbsp;no initial condition<br>
         <strong>else</strong> <strong>der</strong>(p)=0 </td></tr>

<tr><td> SteadyState</td>
    <td> <strong>der</strong>(m)=0  </td>
    <td> no initial conditions </td></tr>
</table>

&nbsp;<br>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><td colspan=\"3\"><strong>Energy balance</strong> </td></tr>
<tr><td><strong>Dynamics.</strong></td>
    <td><strong>Balance equation</strong></td>
    <td><strong>Initial condition</strong></td></tr>

<tr><td> DynamicFreeInitial</td>
    <td> no restrictions </td>
    <td> no initial conditions </td></tr>

<tr><td> FixedInitial</td>
    <td> no restrictions </td>
    <td> T=T_start or h=h_start </td></tr>

<tr><td> SteadyStateInitial</td>
    <td> no restrictions </td>
    <td> <strong>der</strong>(T)=0 or <strong>der</strong>(h)=0 </td></tr>

<tr><td> SteadyState</td>
    <td> <strong>der</strong>(U)=0  </td>
    <td> no initial conditions </td></tr>
</table>

&nbsp;<br>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><td colspan=\"3\"><strong>Momentum balance</strong> </td></tr>
<tr><td><strong>Dynamics.</strong></td>
    <td><strong>Balance equation</strong></td>
    <td><strong>Initial condition</strong></td></tr>

<tr><td> DynamicFreeInitial</td>
    <td> no restrictions </td>
    <td> no initial conditions </td></tr>

<tr><td> FixedInitial</td>
    <td> no restrictions </td>
    <td> m_flow = m_flow_start </td></tr>

<tr><td> SteadyStateInitial</td>
    <td> no restrictions </td>
    <td> <strong>der</strong>(m_flow)=0 </td></tr>

<tr><td> SteadyState</td>
    <td> <strong>der</strong>(m_flow)=0 </td>
    <td> no initial conditions </td></tr>
</table>

<p>
In the tables above, the equations are given for one-substance fluids. For multiple-substance
fluids and for trace substances, equivalent equations hold.
</p>

<p>
Medium.singleState is a medium property and defines whether the medium is only
described by one state (+ the mass fractions in case of a multi-substance fluid). In such
a case one initial condition less must be provided. For example, incompressible
media have Medium.singleState = <strong>true</strong>.
</p>

</html>"));

      type PortFlowDirection = enumeration(
          Entering "Fluid flow is only entering",
          Leaving "Fluid flow is only leaving",
          Bidirectional "No restrictions on fluid flow (flow reversal possible)")
        "Enumeration to define whether flow reversal is allowed" annotation (
          Documentation(info="<html>

<p>
Enumeration to define the assumptions on the model for the
direction of fluid flow at a port (to be selected via choices menu):
</p>

<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
<tr><th><strong>PortFlowDirection.</strong></th>
    <th><strong>Meaning</strong></th></tr>

<tr><td>Entering</td>
    <td>Fluid flow is only entering the port from the outside</td></tr>

<tr><td>Leaving</td>
    <td>Fluid flow is only leaving the port to the outside</td></tr>

<tr><td>Bidirectional</td>
    <td>No restrictions on fluid flow (flow reversal possible)</td></tr>
</table>

<p>
The default is \"PortFlowDirection.Bidirectional\". If you are completely sure that
the flow is only in one direction, then the other settings may
make the simulation of your model faster.
</p>

</html>"));
      annotation (preferredView="info",
                  Documentation(info="<html>

</html>"));
    end Types;

    package Utilities "Utility models to construct fluid components (should not be used directly)"
      extends Modelica.Icons.UtilitiesPackage;

      function checkBoundary "Check whether boundary definition is correct"
        extends Modelica.Icons.Function;
        input String mediumName;
        input String substanceNames[:] "Names of substances";
        input Boolean singleState;
        input Boolean define_p;
        input Real X_boundary[:];
        input String modelName = "??? boundary ???";
    protected
        Integer nX = size(X_boundary,1);
        String X_str;
      algorithm
        assert(not singleState or singleState and define_p, "
Wrong value of parameter define_p (= false) in model \""     + modelName + "\":
The selected medium \""     + mediumName + "\" has Medium.singleState=true.
Therefore, an boundary density cannot be defined and
define_p = true is required.
");

        for i in 1:nX loop
          assert(X_boundary[i] >= 0.0, "
Wrong boundary mass fractions in medium \""
      + mediumName + "\" in model \"" + modelName + "\":
The boundary value X_boundary("   + String(i) + ") = " + String(
            X_boundary[i]) + "
is negative. It must be positive.
");     end for;

        if nX > 0 and abs(sum(X_boundary) - 1.0) > 1e-10 then
           X_str :="";
           for i in 1:nX loop
              X_str :=X_str + "   X_boundary[" + String(i) + "] = " + String(X_boundary[
              i]) + " \"" + substanceNames[i] + "\"\n";
           end for;
           Modelica.Utilities.Streams.error(
              "The boundary mass fractions in medium \"" + mediumName + "\" in model \"" + modelName + "\"\n" +
              "do not sum up to 1. Instead, sum(X_boundary) = " + String(sum(X_boundary)) + ":\n"
              + X_str);
        end if;
      end checkBoundary;

      function regSquare2
        "Anti-symmetric approximation of square with discontinuous factor so that the first derivative is non-zero and is continuous"
        extends Modelica.Icons.Function;
        input Real x "Abscissa value";
        input Real x_small(min=0)=0.01
          "Approximation of function for |x| <= x_small";
        input Real k1(min=0)=1 "y = (if x>=0 then k1 else k2)*x*|x|";
        input Real k2(min=0)=1 "y = (if x>=0 then k1 else k2)*x*|x|";
        input Boolean use_yd0 = false "= true, if yd0 shall be used";
        input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
        output Real y "Ordinate value";
    protected
        encapsulated function regSquare2_utility
          "Interpolating with two 3-order polynomials with a prescribed derivative at x=0"
          import Modelica;
          extends Modelica.Icons.Function;
          import Modelica.Fluid.Utilities.evaluatePoly3_derivativeAtZero;
           input Real x;
           input Real x1 "Approximation of function abs(x) < x1";
           input Real k1 "y = (if x>=0 then k1 else -k2)*x*|x|; k1 >= k2";
           input Real k2 "y = (if x>=0 then k1 else -k2)*x*|x|";
           input Boolean use_yd0 = false "= true, if yd0 shall be used";
           input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
           output Real y;
      protected
           Real x2;
           Real y1;
           Real y2;
           Real y1d;
           Real y2d;
           Real w;
           Real w1;
           Real w2;
           Real y0d;
           Real ww;
        algorithm
           // x2 :=-x1*(k2/k1)^2;
           x2 := -x1;
           if x <= x2 then
              y := -k2*x^2;
           else
               y1 := k1*x1^2;
               y2 :=-k2*x2^2;
              y1d := k1*2*x1;
              y2d :=-k2*2*x2;
              if use_yd0 then
                 y0d :=yd0;
              else
                 /* Determine derivative, such that first and second derivative
              of left and right polynomial are identical at x=0:
              see derivation in function regRoot2
           */
                 w :=x2/x1;
                 y0d := ( (3*y2 - x2*y2d)/w - (3*y1 - x1*y1d)*w) /(2*x1*(1 - w));
              end if;

              /* Modify derivative y0d, such that the polynomial is
           monotonically increasing. A sufficient condition is
             0 <= y0d <= sqrt(5)*k_i*|x_i|
        */
              w1 :=sqrt(5)*k1*x1;
              w2 :=sqrt(5)*k2*abs(x2);
              // y0d :=min(y0d, 0.9*min(w1, w2));
              ww :=0.9*(if w1 < w2 then w1 else w2);
              if ww < y0d then
                 y0d :=ww;
              end if;
              y := if x >= 0 then evaluatePoly3_derivativeAtZero(x,x1,y1,y1d,y0d) else
                                  evaluatePoly3_derivativeAtZero(x,x2,y2,y2d,y0d);
           end if;
           annotation(smoothOrder=2);
        end regSquare2_utility;
      algorithm
        y := smooth(2,if x >= x_small then k1*x^2 else
                      if x <= -x_small then -k2*x^2 else
                      if k1 >= k2 then regSquare2_utility(x,x_small,k1,k2,use_yd0,yd0) else
                                      -regSquare2_utility(-x,x_small,k2,k1,use_yd0,yd0));
        annotation(smoothOrder=2, Documentation(info="<html>
<p>
Approximates the function
</p>
<blockquote><pre>
y = <strong>if</strong> x &ge; 0 <strong>then</strong> k1*x*x <strong>else</strong> -k2*x*x, with k1, k2 > 0
</pre></blockquote>
<p>
in such a way that within the region -x_small &le; x &le; x_small,
the function is described by two polynomials of third order
(one in the region -x_small .. 0 and one within the region 0 .. x_small)
such that
</p>

<ul>
<li> The derivative at x=0 is non-zero (in order that the
     inverse of the function does not have an infinite derivative).</li>
<li> The overall function is continuous with a
     continuous first derivative everywhere.</li>
<li> If parameter use_yd0 = <strong>false</strong>, the two polynomials
     are constructed such that the second derivatives at x=0
     are identical. If use_yd0 = <strong>true</strong>, the derivative
     at x=0 is explicitly provided via the additional argument
     yd0. If necessary, the derivative yd0 is automatically
     reduced in order that the polynomials are strict monotonically
     increasing <em>[Fritsch and Carlson, 1980]</em>.</li>
</ul>

<p>
A typical screenshot for k1=1, k2=3 is shown in the next figure:
</p>
<p>
<img src=\"modelica://Modelica/Resources/Images/Fluid/Utilities/regSquare2_b.png\"
     alt=\"regSquare2_b.png\">
</p>

<p>
The (smooth, non-zero) derivative of the function with
k1=1, k2=3 is shown in the next figure:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Fluid/Utilities/regSquare2_c.png\"
     alt=\"regSquare2_b.png\">
</p>

<p>
<strong>Literature</strong>
</p>

<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <strong>Monotone piecewise cubic interpolation</strong>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>",     revisions="<html>
<ul>
<li><em>Nov., 2005</em>
    by <a href=\"mailto:Martin.Otter@DLR.de\">Martin Otter</a>:<br>
    Designed and implemented.</li>
</ul>
</html>"));
      end regSquare2;

      function regStep
        "Approximation of a general step, such that the characteristic is continuous and differentiable"
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
        annotation(Documentation(revisions="<html>
<ul>
<li><em>April 29, 2008</em>
    by <a href=\"mailto:Martin.Otter@DLR.de\">Martin Otter</a>:<br>
    Designed and implemented.</li>
<li><em>August 12, 2008</em>
    by <a href=\"mailto:Michael.Sielemann@dlr.de\">Michael Sielemann</a>:<br>
    Minor modification to cover the limit case <code>x_small -> 0</code> without division by zero.</li>
</ul>
</html>",       info="<html>
<p>
This function is used to approximate the equation
</p>
<blockquote><pre>
y = <strong>if</strong> x &gt; 0 <strong>then</strong> y1 <strong>else</strong> y2;
</pre></blockquote>

<p>
by a smooth characteristic, so that the expression is continuous and differentiable:
</p>

<blockquote><pre>
y = <strong>smooth</strong>(1, <strong>if</strong> x &gt;  x_small <strong>then</strong> y1 <strong>else</strong>
              <strong>if</strong> x &lt; -x_small <strong>then</strong> y2 <strong>else</strong> f(y1, y2));
</pre></blockquote>

<p>
In the region -x_small &lt; x &lt; x_small a 2nd order polynomial is used
for a smooth transition from y1 to y2.
</p>
</html>"));
      end regStep;

      function evaluatePoly3_derivativeAtZero
        "Evaluate polynomial of order 3 that passes the origin with a predefined derivative"
        extends Modelica.Icons.Function;
        input Real x "Value for which polynomial shall be evaluated";
        input Real x1 "Abscissa value";
        input Real y1 "y1=f(x1)";
        input Real y1d "First derivative at y1";
        input Real y0d "First derivative at f(x=0)";
        output Real y;
    protected
        Real a1;
        Real a2;
        Real a3;
        Real xx;
      algorithm
        a1 := x1*y0d;
        a2 := 3*y1 - x1*y1d - 2*a1;
        a3 := y1 - a2 - a1;
        xx := x/x1;
        y  := xx*(a1 + xx*(a2 + xx*a3));
        annotation(smoothOrder=3, Documentation(info="<html>

</html>"));
      end evaluatePoly3_derivativeAtZero;
      annotation (Documentation(info="<html>

</html>"));
    end Utilities;
  annotation (Icon(graphics={
          Polygon(points={{-70,26},{68,-44},{68,26},{2,-10},{-70,-42},{-70,26}}),
          Line(points={{2,42},{2,-10}}),
          Rectangle(
            extent={{-18,50},{22,42}},
            fillPattern=FillPattern.Solid)}), preferredView="info",
    Documentation(info="<html>
<p>
Library <strong>Modelica.Fluid</strong> is a <strong>free</strong> Modelica package providing components for
<strong>1-dimensional thermo-fluid flow</strong> in networks of vessels, pipes, fluid machines, valves and fittings.
A unique feature is that the component equations and the media models
as well as pressure loss and heat transfer correlations are decoupled from each other.
All components are implemented such that they can be used for
media from the Modelica.Media library. This means especially that an
incompressible or compressible medium, a single or a multiple
substance medium with one or more phases might be used.
</p>

<p>
In the next figure, several features of the library are demonstrated with
a simple heating system with a closed flow cycle. By just changing one configuration parameter in the system object the equations are changed between steady-state and dynamic simulation with fixed or steady-state initial conditions.
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Fluid/HeatingSystem.png\" border=\"1\"
     alt=\"HeatingSystem.png\">
</p>

<p>
With respect to previous versions, the design
of the connectors has been changed in a non-backward compatible way,
using the recently developed concept
of stream connectors that results in much more reliable simulations
(see also <a href=\"modelica://Modelica/Resources/Documentation/Fluid/Stream-Connectors-Overview-Rationale.pdf\">Stream-Connectors-Overview-Rationale.pdf</a>).
This extension was included in Modelica 3.1.
</p>

<p>
The following parts are useful, when newly starting with this library:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Fluid.UsersGuide\">Modelica.Fluid.UsersGuide</a>.</li>
<li> <a href=\"modelica://Modelica.Fluid.UsersGuide.ReleaseNotes\">Modelica.Fluid.UsersGuide.ReleaseNotes</a>
     summarizes the changes of the library releases.</li>
<li> <a href=\"modelica://Modelica.Fluid.Examples\">Modelica.Fluid.Examples</a>
     contains examples that demonstrate the usage of this library.</li>
</ul>
<p>
Copyright &copy; 2002-2020, Modelica Association and contributors
</p>
</html>"));
  end Fluid;

  package Media "Library of media property models"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;
    import Cv = Modelica.Units.Conversions;

  package Interfaces "Interfaces for media models"
    extends Modelica.Icons.InterfacesPackage;

    partial package PartialMedium
      "Partial medium properties (base package of all media packages)"
      extends Modelica.Media.Interfaces.Types;
      extends Modelica.Icons.MaterialPropertiesPackage;

      // Constants to be set in Medium
      constant Modelica.Media.Interfaces.Choices.IndependentVariables
        ThermoStates "Enumeration type for independent variables";
      constant String mediumName="unusablePartialMedium" "Name of the medium";
      constant String substanceNames[:]={mediumName}
        "Names of the mixture substances. Set substanceNames={mediumName} if only one substance.";
      constant String extraPropertiesNames[:]=fill("", 0)
        "Names of the additional (extra) transported properties. Set extraPropertiesNames=fill(\"\",0) if unused";
      constant Boolean singleState
        "= true, if u and d are not a function of pressure";
      constant Boolean reducedX=true
        "= true if medium contains the equation sum(X) = 1.0; set reducedX=true if only one substance (see docu for details)";
      constant Boolean fixedX=false
        "= true if medium contains the equation X = reference_X";
      constant AbsolutePressure reference_p=101325
        "Reference pressure of Medium: default 1 atmosphere";
      constant Temperature reference_T=298.15
        "Reference temperature of Medium: default 25 deg Celsius";
      constant MassFraction reference_X[nX]=fill(1/nX, nX)
        "Default mass fractions of medium";
      constant AbsolutePressure p_default=101325
        "Default value for pressure of medium (for initialization)";
      constant Temperature T_default=Modelica.Units.Conversions.from_degC(20)
        "Default value for temperature of medium (for initialization)";
      constant SpecificEnthalpy h_default=specificEnthalpy_pTX(
              p_default,
              T_default,
              X_default)
        "Default value for specific enthalpy of medium (for initialization)";
      constant MassFraction X_default[nX]=reference_X
        "Default value for mass fractions of medium (for initialization)";
      constant ExtraProperty C_default[nC]=fill(0, nC)
        "Default value for trace substances of medium (for initialization)";

      final constant Integer nS=size(substanceNames, 1) "Number of substances";
      constant Integer nX=nS "Number of mass fractions";
      constant Integer nXi=if fixedX then 0 else if reducedX then nS - 1 else nS
        "Number of structurally independent mass fractions (see docu for details)";

      final constant Integer nC=size(extraPropertiesNames, 1)
        "Number of extra (outside of standard mass-balance) transported properties";
      constant Real C_nominal[nC](min=fill(Modelica.Constants.eps, nC)) = 1.0e-6*
        ones(nC) "Default for the nominal values for the extra properties";
      replaceable record FluidConstants =
          Modelica.Media.Interfaces.Types.Basic.FluidConstants
        "Critical, triple, molecular and other standard data of fluid";

      replaceable record ThermodynamicState
        "Minimal variable set that is available as input argument to every medium function"
        extends Modelica.Icons.Record;
      end ThermodynamicState;

      replaceable partial model BaseProperties
        "Base properties (p, d, T, h, u, R_s, MM and, if applicable, X and Xi) of a medium"
        InputAbsolutePressure p "Absolute pressure of medium";
        InputMassFraction[nXi] Xi(start=reference_X[1:nXi])
          "Structurally independent mass fractions";
        InputSpecificEnthalpy h "Specific enthalpy of medium";
        Density d "Density of medium";
        Temperature T "Temperature of medium";
        MassFraction[nX] X(start=reference_X)
          "Mass fractions (= (component mass)/total mass  m_i/m)";
        SpecificInternalEnergy u "Specific internal energy of medium";
        SpecificHeatCapacity R_s "Gas constant (of mixture if applicable)";
        MolarMass MM "Molar mass (of mixture or single fluid)";
        ThermodynamicState state
          "Thermodynamic state record for optional functions";
        parameter Boolean preferredMediumStates=false
          "= true if StateSelect.prefer shall be used for the independent property variables of the medium"
          annotation (Evaluate=true, Dialog(tab="Advanced"));
        parameter Boolean standardOrderComponents=true
          "If true, and reducedX = true, the last element of X will be computed from the other ones";
        Modelica.Units.NonSI.Temperature_degC T_degC=
            Modelica.Units.Conversions.to_degC(T)
          "Temperature of medium in [degC]";
        Modelica.Units.NonSI.Pressure_bar p_bar=
            Modelica.Units.Conversions.to_bar(p)
          "Absolute pressure of medium in [bar]";

        // Local connector definition, used for equation balancing check
        connector InputAbsolutePressure = input SI.AbsolutePressure
          "Pressure as input signal connector";
        connector InputSpecificEnthalpy = input SI.SpecificEnthalpy
          "Specific enthalpy as input signal connector";
        connector InputMassFraction = input SI.MassFraction
          "Mass fraction as input signal connector";

      equation
        if standardOrderComponents then
          Xi = X[1:nXi];

          if fixedX then
            X = reference_X;
          end if;
          if reducedX and not fixedX then
            X[nX] = 1 - sum(Xi);
          end if;
          for i in 1:nX loop
            assert(X[i] >= -1.e-5 and X[i] <= 1 + 1.e-5, "Mass fraction X[" +
              String(i) + "] = " + String(X[i]) + "of substance " +
              substanceNames[i] + "\nof medium " + mediumName +
              " is not in the range 0..1");
          end for;

        end if;

        assert(p >= 0.0, "Pressure (= " + String(p) + " Pa) of medium \"" +
          mediumName + "\" is negative\n(Temperature = " + String(T) + " K)");
        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Rectangle(
                extent={{-100,100},{100,-100}},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid,
                lineColor={0,0,255}), Text(
                extent={{-152,164},{152,102}},
                textString="%name",
                textColor={0,0,255})}), Documentation(info="<html>
<p>
Model <strong>BaseProperties</strong> is a model within package <strong>PartialMedium</strong>
and contains the <strong>declarations</strong> of the minimum number of
variables that every medium model is supposed to support.
A specific medium inherits from model <strong>BaseProperties</strong> and provides
the equations for the basic properties.</p>
<p>
The BaseProperties model contains the following <strong>7+nXi variables</strong>
(nXi is the number of independent mass fractions defined in package
PartialMedium):
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>Variable</strong></td>
      <td><strong>Unit</strong></td>
      <td><strong>Description</strong></td></tr>
  <tr><td>T</td>
      <td>K</td>
      <td>Temperature</td></tr>
  <tr><td>p</td>
      <td>Pa</td>
      <td>Absolute pressure</td></tr>
  <tr><td>d</td>
      <td>kg/m3</td>
      <td>Density</td></tr>
  <tr><td>h</td>
      <td>J/kg</td>
      <td>Specific enthalpy</td></tr>
  <tr><td>u</td>
      <td>J/kg</td>
      <td>Specific internal energy</td></tr>
  <tr><td>Xi[nXi]</td>
      <td>kg/kg</td>
      <td>Structurally independent mass fractions</td></tr>
  <tr><td>R_s</td>
      <td>J/(kg.K)</td>
      <td>Specific gas constant (of mixture if applicable)</td></tr>
  <tr><td>MM</td>
      <td>kg/mol</td>
      <td>Molar mass</td></tr>
</table>
<p>
In order to implement an actual medium model, one can extend from this
base model and add <strong>5 equations</strong> that provide relations among
these variables. Equations will also have to be added in order to
set all the variables within the ThermodynamicState record state.</p>
<p>
If standardOrderComponents=true, the full composition vector X[nX]
is determined by the equations contained in this base class, depending
on the independent mass fraction vector Xi[nXi].</p>
<p>Additional <strong>2 + nXi</strong> equations will have to be provided
when using the BaseProperties model, in order to fully specify the
thermodynamic conditions. The input connector qualifier applied to
p, h, and nXi indirectly declares the number of missing equations,
permitting advanced equation balance checking by Modelica tools.
Please note that this doesn't mean that the additional equations
should be connection equations, nor that exactly those variables
should be supplied, in order to complete the model.
For further information, see the <a href=\"modelica://Modelica.Media.UsersGuide\">Modelica.Media User's guide</a>, and
<a href=\"https://specification.modelica.org/v3.4/Ch4.html#balanced-models\">Section 4.7 (Balanced Models) of the Modelica 3.4 specification</a>.</p>
</html>"));
      end BaseProperties;

      replaceable partial function setState_pTX
        "Return thermodynamic state as function of p, T and composition X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_pTX;

      replaceable partial function setState_phX
        "Return thermodynamic state as function of p, h and composition X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_phX;

      replaceable partial function setState_psX
        "Return thermodynamic state as function of p, s and composition X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_psX;

      replaceable partial function setState_dTX
        "Return thermodynamic state as function of d, T and composition X or Xi"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      end setState_dTX;

      replaceable partial function setSmoothState
        "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
        extends Modelica.Icons.Function;
        input Real x "m_flow or dp";
        input ThermodynamicState state_a "Thermodynamic state if x > 0";
        input ThermodynamicState state_b "Thermodynamic state if x < 0";
        input Real x_small(min=0)
          "Smooth transition in the region -x_small < x < x_small";
        output ThermodynamicState state
          "Smooth thermodynamic state for all x (continuous and differentiable)";
        annotation (Documentation(info="<html>
<p>
This function is used to approximate the equation
</p>
<blockquote><pre>
state = <strong>if</strong> x &gt; 0 <strong>then</strong> state_a <strong>else</strong> state_b;
</pre></blockquote>

<p>
by a smooth characteristic, so that the expression is continuous and differentiable:
</p>

<blockquote><pre>
state := <strong>smooth</strong>(1, <strong>if</strong> x &gt;  x_small <strong>then</strong> state_a <strong>else</strong>
                   <strong>if</strong> x &lt; -x_small <strong>then</strong> state_b <strong>else</strong> f(state_a, state_b));
</pre></blockquote>

<p>
This is performed by applying function <strong>Media.Common.smoothStep</strong>(..)
on every element of the thermodynamic state record.
</p>

<p>
If <strong>mass fractions</strong> X[:] are approximated with this function then this can be performed
for all <strong>nX</strong> mass fractions, instead of applying it for nX-1 mass fractions and computing
the last one by the mass fraction constraint sum(X)=1. The reason is that the approximating function has the
property that sum(state.X) = 1, provided sum(state_a.X) = sum(state_b.X) = 1.
This can be shown by evaluating the approximating function in the abs(x) &lt; x_small
region (otherwise state.X is either state_a.X or state_b.X):
</p>

<blockquote><pre>
X[1]  = smoothStep(x, X_a[1] , X_b[1] , x_small);
X[2]  = smoothStep(x, X_a[2] , X_b[2] , x_small);
   ...
X[nX] = smoothStep(x, X_a[nX], X_b[nX], x_small);
</pre></blockquote>

<p>
or
</p>

<blockquote><pre>
X[1]  = c*(X_a[1]  - X_b[1])  + (X_a[1]  + X_b[1])/2
X[2]  = c*(X_a[2]  - X_b[2])  + (X_a[2]  + X_b[2])/2;
   ...
X[nX] = c*(X_a[nX] - X_b[nX]) + (X_a[nX] + X_b[nX])/2;
c     = (x/x_small)*((x/x_small)^2 - 3)/4
</pre></blockquote>

<p>
Summing all mass fractions together results in
</p>

<blockquote><pre>
sum(X) = c*(sum(X_a) - sum(X_b)) + (sum(X_a) + sum(X_b))/2
       = c*(1 - 1) + (1 + 1)/2
       = 1
</pre></blockquote>

</html>"));
      end setSmoothState;

      replaceable partial function dynamicViscosity "Return dynamic viscosity"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DynamicViscosity eta "Dynamic viscosity";
      end dynamicViscosity;

      replaceable partial function thermalConductivity
        "Return thermal conductivity"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output ThermalConductivity lambda "Thermal conductivity";
      end thermalConductivity;

      replaceable function prandtlNumber "Return the Prandtl number"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output PrandtlNumber Pr "Prandtl number";
      algorithm
        Pr := dynamicViscosity(state)*specificHeatCapacityCp(state)/
          thermalConductivity(state);
      end prandtlNumber;

      replaceable partial function pressure "Return pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output AbsolutePressure p "Pressure";
      end pressure;

      replaceable partial function temperature "Return temperature"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Temperature T "Temperature";
      end temperature;

      replaceable partial function density "Return density"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Density d "Density";
      end density;

      replaceable partial function specificEnthalpy "Return specific enthalpy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnthalpy h "Specific enthalpy";
      end specificEnthalpy;

      replaceable partial function specificInternalEnergy
        "Return specific internal energy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnergy u "Specific internal energy";
      end specificInternalEnergy;

      replaceable partial function specificEntropy "Return specific entropy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEntropy s "Specific entropy";
      end specificEntropy;

      replaceable partial function specificGibbsEnergy
        "Return specific Gibbs energy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnergy g "Specific Gibbs energy";
      end specificGibbsEnergy;

      replaceable partial function specificHelmholtzEnergy
        "Return specific Helmholtz energy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificEnergy f "Specific Helmholtz energy";
      end specificHelmholtzEnergy;

      replaceable partial function specificHeatCapacityCp
        "Return specific heat capacity at constant pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificHeatCapacity cp
          "Specific heat capacity at constant pressure";
      end specificHeatCapacityCp;

      function heatCapacity_cp = specificHeatCapacityCp
        "Alias for deprecated name";

      replaceable partial function specificHeatCapacityCv
        "Return specific heat capacity at constant volume"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SpecificHeatCapacity cv
          "Specific heat capacity at constant volume";
      end specificHeatCapacityCv;

      function heatCapacity_cv = specificHeatCapacityCv
        "Alias for deprecated name";

      replaceable partial function isentropicExponent
        "Return isentropic exponent"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output IsentropicExponent gamma "Isentropic exponent";
      end isentropicExponent;

      replaceable partial function isentropicEnthalpy
        "Return isentropic enthalpy"
        extends Modelica.Icons.Function;
        input AbsolutePressure p_downstream "Downstream pressure";
        input ThermodynamicState refState "Reference state for entropy";
        output SpecificEnthalpy h_is "Isentropic enthalpy";
        annotation (Documentation(info="<html>
<p>
This function computes an isentropic state transformation:
</p>
<ol>
<li> A medium is in a particular state, refState.</li>
<li> The enthalpy at another state (h_is) shall be computed
     under the assumption that the state transformation from refState to h_is
     is performed with a change of specific entropy ds = 0 and the pressure of state h_is
     is p_downstream and the composition X upstream and downstream is assumed to be the same.</li>
</ol>

</html>"));
      end isentropicEnthalpy;

      replaceable partial function velocityOfSound "Return velocity of sound"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output VelocityOfSound a "Velocity of sound";
      end velocityOfSound;

      replaceable partial function isobaricExpansionCoefficient
        "Return overall the isobaric expansion coefficient beta"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output IsobaricExpansionCoefficient beta "Isobaric expansion coefficient";
        annotation (Documentation(info="<html>
<blockquote><pre>
beta is defined as  1/v * der(v,T), with v = 1/d, at constant pressure p.
</pre></blockquote>
</html>"));
      end isobaricExpansionCoefficient;

      function beta = isobaricExpansionCoefficient
        "Alias for isobaricExpansionCoefficient for user convenience";

      replaceable partial function isothermalCompressibility
        "Return overall the isothermal compressibility factor"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output SI.IsothermalCompressibility kappa "Isothermal compressibility";
        annotation (Documentation(info="<html>
<blockquote><pre>

kappa is defined as - 1/v * der(v,p), with v = 1/d at constant temperature T.

</pre></blockquote>
</html>"));
      end isothermalCompressibility;

      function kappa = isothermalCompressibility
        "Alias of isothermalCompressibility for user convenience";

      // explicit derivative functions for finite element models
      replaceable partial function density_derp_h
        "Return density derivative w.r.t. pressure at const specific enthalpy"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByPressure ddph "Density derivative w.r.t. pressure";
      end density_derp_h;

      replaceable partial function density_derh_p
        "Return density derivative w.r.t. specific enthalpy at constant pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByEnthalpy ddhp
          "Density derivative w.r.t. specific enthalpy";
      end density_derh_p;

      replaceable partial function density_derp_T
        "Return density derivative w.r.t. pressure at const temperature"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByPressure ddpT "Density derivative w.r.t. pressure";
      end density_derp_T;

      replaceable partial function density_derT_p
        "Return density derivative w.r.t. temperature at constant pressure"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output DerDensityByTemperature ddTp
          "Density derivative w.r.t. temperature";
      end density_derT_p;

      replaceable partial function density_derX
        "Return density derivative w.r.t. mass fraction"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output Density[nX] dddX "Derivative of density w.r.t. mass fraction";
      end density_derX;

      replaceable partial function molarMass
        "Return the molar mass of the medium"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state record";
        output MolarMass MM "Mixture molar mass";
      end molarMass;

      replaceable function specificEnthalpy_pTX
        "Return specific enthalpy from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_pTX(
                p,
                T,
                X));
        annotation (inverse(T=temperature_phX(
                      p,
                      h,
                      X)));
      end specificEnthalpy_pTX;

      replaceable function specificEntropy_pTX
        "Return specific enthalpy from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output SpecificEntropy s "Specific entropy";
      algorithm
        s := specificEntropy(setState_pTX(
                p,
                T,
                X));

        annotation (inverse(T=temperature_psX(
                      p,
                      s,
                      X)));
      end specificEntropy_pTX;

      replaceable function density_pTX "Return density from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:] "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_pTX(
                p,
                T,
                X));
      end density_pTX;

      replaceable function temperature_phX
        "Return temperature from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_phX(
                p,
                h,
                X));
      end temperature_phX;

      replaceable function density_phX "Return density from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_phX(
                p,
                h,
                X));
      end density_phX;

      replaceable function temperature_psX
        "Return temperature from p,s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Temperature T "Temperature";
      algorithm
        T := temperature(setState_psX(
                p,
                s,
                X));
        annotation (inverse(s=specificEntropy_pTX(
                      p,
                      T,
                      X)));
      end temperature_psX;

      replaceable function density_psX "Return density from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_psX(
                p,
                s,
                X));
      end density_psX;

      replaceable function specificEnthalpy_psX
        "Return specific enthalpy from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_psX(
                p,
                s,
                X));
      end specificEnthalpy_psX;

      type MassFlowRate = SI.MassFlowRate (
          quantity="MassFlowRate." + mediumName,
          min=-1.0e5,
          max=1.e5) "Type for mass flow rate with medium specific attributes";

      annotation (Documentation(info="<html>
<p>
<strong>PartialMedium</strong> is a package and contains all <strong>declarations</strong> for
a medium. This means that constants, models, and functions
are defined that every medium is supposed to support
(some of them are optional). A medium package
inherits from <strong>PartialMedium</strong> and provides the
equations for the medium. The details of this package
are described in
<a href=\"modelica://Modelica.Media.UsersGuide\">Modelica.Media.UsersGuide</a>.
</p>
</html>",   revisions="<html>

</html>"));
    end PartialMedium;

    partial package PartialPureSubstance
      "Base class for pure substances of one chemical substance"
      extends PartialMedium(final reducedX=true, final fixedX=true);

      replaceable function setState_pT "Return thermodynamic state from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_pTX(
                p,
                T,
                fill(0, 0));
      end setState_pT;

      replaceable function setState_ph "Return thermodynamic state from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_phX(
                p,
                h,
                fill(0, 0));
      end setState_ph;

      replaceable function setState_ps "Return thermodynamic state from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_psX(
                p,
                s,
                fill(0, 0));
      end setState_ps;

      replaceable function setState_dT "Return thermodynamic state from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := setState_dTX(
                d,
                T,
                fill(0, 0));
      end setState_dT;

      replaceable function density_ph "Return density from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        output Density d "Density";
      algorithm
        d := density_phX(
                p,
                h,
                fill(0, 0));
      end density_ph;

      replaceable function temperature_ph "Return temperature from p and h"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        output Temperature T "Temperature";
      algorithm
        T := temperature_phX(
                p,
                h,
                fill(0, 0));
      end temperature_ph;

      replaceable function pressure_dT "Return pressure from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        output AbsolutePressure p "Pressure";
      algorithm
        p := pressure(setState_dTX(
                d,
                T,
                fill(0, 0)));
      end pressure_dT;

      replaceable function specificEnthalpy_dT
        "Return specific enthalpy from d and T"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy(setState_dTX(
                d,
                T,
                fill(0, 0)));
      end specificEnthalpy_dT;

      replaceable function specificEnthalpy_ps
        "Return specific enthalpy from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy_psX(
                p,
                s,
                fill(0, 0));
      end specificEnthalpy_ps;

      replaceable function temperature_ps "Return temperature from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        output Temperature T "Temperature";
      algorithm
        T := temperature_psX(
                p,
                s,
                fill(0, 0));
      end temperature_ps;

      replaceable function density_ps "Return density from p and s"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        output Density d "Density";
      algorithm
        d := density_psX(
                p,
                s,
                fill(0, 0));
      end density_ps;

      replaceable function specificEnthalpy_pT
        "Return specific enthalpy from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := specificEnthalpy_pTX(
                p,
                T,
                fill(0, 0));
      end specificEnthalpy_pT;

      replaceable function density_pT "Return density from p and T"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        output Density d "Density";
      algorithm
        d := density(setState_pTX(
                p,
                T,
                fill(0, 0)));
      end density_pT;

      redeclare replaceable partial model extends BaseProperties(final
          standardOrderComponents=true)
      end BaseProperties;
    end PartialPureSubstance;

    partial package PartialMixtureMedium
      "Base class for pure substances of several chemical substances"
      extends PartialMedium(redeclare replaceable record FluidConstants =
            Modelica.Media.Interfaces.Types.IdealGas.FluidConstants);

      redeclare replaceable record extends ThermodynamicState
        "Thermodynamic state variables"
        AbsolutePressure p "Absolute pressure of medium";
        Temperature T "Temperature of medium";
        MassFraction[nX] X(start=reference_X)
          "Mass fractions (= (component mass)/total mass  m_i/m)";
      end ThermodynamicState;

      constant FluidConstants[nS] fluidConstants "Constant data for the fluid";

      replaceable function gasConstant
        "Return the gas constant of the mixture (also for liquids)"
        extends Modelica.Icons.Function;
        input ThermodynamicState state "Thermodynamic state";
        output SI.SpecificHeatCapacity R_s "Mixture gas constant";
      end gasConstant;

      function moleToMassFractions "Return mass fractions X from mole fractions"
        extends Modelica.Icons.Function;
        input SI.MoleFraction moleFractions[:] "Mole fractions of mixture";
        input MolarMass[:] MMX "Molar masses of components";
        output SI.MassFraction X[size(moleFractions, 1)]
          "Mass fractions of gas mixture";
    protected
        MolarMass Mmix=moleFractions*MMX "Molar mass of mixture";
      algorithm
        for i in 1:size(moleFractions, 1) loop
          X[i] := moleFractions[i]*MMX[i]/Mmix;
        end for;
        annotation (smoothOrder=5);
      end moleToMassFractions;

      function massToMoleFractions "Return mole fractions from mass fractions X"
        extends Modelica.Icons.Function;
        input SI.MassFraction X[:] "Mass fractions of mixture";
        input SI.MolarMass[:] MMX "Molar masses of components";
        output SI.MoleFraction moleFractions[size(X, 1)]
          "Mole fractions of gas mixture";
    protected
        Real invMMX[size(X, 1)] "Inverses of molar weights";
        SI.MolarMass Mmix "Molar mass of mixture";
      algorithm
        for i in 1:size(X, 1) loop
          invMMX[i] := 1/MMX[i];
        end for;
        Mmix := 1/(X*invMMX);
        for i in 1:size(X, 1) loop
          moleFractions[i] := Mmix*X[i]/MMX[i];
        end for;
        annotation (smoothOrder=5);
      end massToMoleFractions;

    end PartialMixtureMedium;

    partial package PartialCondensingGases
      "Base class for mixtures of condensing and non-condensing gases"
      extends PartialMixtureMedium(ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pTX);

      replaceable partial function saturationPressure
        "Return saturation pressure of condensing fluid"
        extends Modelica.Icons.Function;
        input Temperature Tsat "Saturation temperature";
        output AbsolutePressure psat "Saturation pressure";
      end saturationPressure;

      replaceable partial function enthalpyOfVaporization
        "Return vaporization enthalpy of condensing fluid"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        output SpecificEnthalpy r0 "Vaporization enthalpy";
      end enthalpyOfVaporization;

      replaceable partial function enthalpyOfLiquid
        "Return liquid enthalpy of condensing fluid"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        output SpecificEnthalpy h "Liquid enthalpy";
      end enthalpyOfLiquid;

      replaceable partial function enthalpyOfGas
        "Return enthalpy of non-condensing gas mixture"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        input MassFraction[:] X "Vector of mass fractions";
        output SpecificEnthalpy h "Specific enthalpy";
      end enthalpyOfGas;

      replaceable partial function enthalpyOfCondensingGas
        "Return enthalpy of condensing gas (most often steam)"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        output SpecificEnthalpy h "Specific enthalpy";
      end enthalpyOfCondensingGas;

      replaceable partial function enthalpyOfNonCondensingGas
        "Return enthalpy of the non-condensing species"
        extends Modelica.Icons.Function;
        input Temperature T "Temperature";
        output SpecificEnthalpy h "Specific enthalpy";
      end enthalpyOfNonCondensingGas;
    end PartialCondensingGases;

    partial package PartialSimpleMedium
      "Medium model with linear dependency of u, h from temperature. All other quantities, especially density, are constant."

      extends Interfaces.PartialPureSubstance(final ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pT,
          final singleState=true);

      constant SpecificHeatCapacity cp_const
        "Constant specific heat capacity at constant pressure";
      constant SpecificHeatCapacity cv_const
        "Constant specific heat capacity at constant volume";
      constant Density d_const "Constant density";
      constant DynamicViscosity eta_const "Constant dynamic viscosity";
      constant ThermalConductivity lambda_const "Constant thermal conductivity";
      constant VelocityOfSound a_const "Constant velocity of sound";
      constant Temperature T_min "Minimum temperature valid for medium model";
      constant Temperature T_max "Maximum temperature valid for medium model";
      constant Temperature T0=reference_T "Zero enthalpy temperature";
      constant MolarMass MM_const "Molar mass";

      constant FluidConstants[nS] fluidConstants "Fluid constants";

      redeclare record extends ThermodynamicState "Thermodynamic state"
        AbsolutePressure p "Absolute pressure of medium";
        Temperature T "Temperature of medium";
      end ThermodynamicState;

      redeclare replaceable model extends BaseProperties(T(stateSelect=if
              preferredMediumStates then StateSelect.prefer else StateSelect.default),
          p(stateSelect=if preferredMediumStates then StateSelect.prefer else
              StateSelect.default)) "Base properties"
      equation
        assert(T >= T_min and T <= T_max, "
Temperature T (= "   + String(T) + " K) is not
in the allowed range ("   + String(T_min) + " K <= T <= " + String(T_max) + " K)
required from medium model \""   + mediumName + "\".
");

        // h = cp_const*(T-T0);
        h = specificEnthalpy_pTX(
                p,
                T,
                X);
        u = cv_const*(T - T0);
        d = d_const;
        R_s = 0;
        MM = MM_const;
        state.T = T;
        state.p = p;
        annotation (Documentation(info="<html>
<p>
This is the most simple incompressible medium model, where
specific enthalpy h and specific internal energy u are only
a function of temperature T and all other provided medium
quantities are assumed to be constant.
Note that the (small) influence of the pressure term p/d is neglected.
</p>
</html>"));
      end BaseProperties;

      redeclare function setState_pTX
        "Return thermodynamic state from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := ThermodynamicState(p=p, T=T);
      end setState_pTX;

      redeclare function setState_phX
        "Return thermodynamic state from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := ThermodynamicState(p=p, T=T0 + h/cp_const);
      end setState_phX;

      redeclare replaceable function setState_psX
        "Return thermodynamic state from p, s, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEntropy s "Specific entropy";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        state := ThermodynamicState(p=p, T=Modelica.Math.exp(s/cp_const +
          Modelica.Math.log(reference_T)))
          "Here the incompressible limit is used, with cp as heat capacity";
      end setState_psX;

      redeclare function setState_dTX
        "Return thermodynamic state from d, T, and X or Xi"
        extends Modelica.Icons.Function;
        input Density d "Density";
        input Temperature T "Temperature";
        input MassFraction X[:]=reference_X "Mass fractions";
        output ThermodynamicState state "Thermodynamic state record";
      algorithm
        assert(false,
          "Pressure can not be computed from temperature and density for an incompressible fluid!");
      end setState_dTX;

      redeclare function extends setSmoothState
        "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
      algorithm
        state := ThermodynamicState(p=Media.Common.smoothStep(
                x,
                state_a.p,
                state_b.p,
                x_small), T=Media.Common.smoothStep(
                x,
                state_a.T,
                state_b.T,
                x_small));
      end setSmoothState;

      redeclare function extends dynamicViscosity "Return dynamic viscosity"

      algorithm
        eta := eta_const;
      end dynamicViscosity;

      redeclare function extends thermalConductivity
        "Return thermal conductivity"

      algorithm
        lambda := lambda_const;
      end thermalConductivity;

      redeclare function extends pressure "Return pressure"

      algorithm
        p := state.p;
      end pressure;

      redeclare function extends temperature "Return temperature"

      algorithm
        T := state.T;
      end temperature;

      redeclare function extends density "Return density"

      algorithm
        d := d_const;
      end density;

      redeclare function extends specificEnthalpy "Return specific enthalpy"

      algorithm
        h := cp_const*(state.T - T0);
      end specificEnthalpy;

      redeclare function extends specificHeatCapacityCp
        "Return specific heat capacity at constant pressure"

      algorithm
        cp := cp_const;
      end specificHeatCapacityCp;

      redeclare function extends specificHeatCapacityCv
        "Return specific heat capacity at constant volume"

      algorithm
        cv := cv_const;
      end specificHeatCapacityCv;

      redeclare function extends isentropicExponent "Return isentropic exponent"

      algorithm
        gamma := cp_const/cv_const;
      end isentropicExponent;

      redeclare function extends velocityOfSound "Return velocity of sound"

      algorithm
        a := a_const;
      end velocityOfSound;

      redeclare function specificEnthalpy_pTX
        "Return specific enthalpy from p, T, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input Temperature T "Temperature";
        input MassFraction X[nX] "Mass fractions";
        output SpecificEnthalpy h "Specific enthalpy";
      algorithm
        h := cp_const*(T - T0);
        annotation (Documentation(info="<html>
<p>
This function computes the specific enthalpy of the fluid, but neglects the (small) influence of the pressure term p/d.
</p>
</html>"));
      end specificEnthalpy_pTX;

      redeclare function temperature_phX
        "Return temperature from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[nX] "Mass fractions";
        output Temperature T "Temperature";
      algorithm
        T := T0 + h/cp_const;
      end temperature_phX;

      redeclare function density_phX "Return density from p, h, and X or Xi"
        extends Modelica.Icons.Function;
        input AbsolutePressure p "Pressure";
        input SpecificEnthalpy h "Specific enthalpy";
        input MassFraction X[nX] "Mass fractions";
        output Density d "Density";
      algorithm
        d := density(setState_phX(
                p,
                h,
                X));
      end density_phX;

      redeclare function extends specificInternalEnergy
        "Return specific internal energy"
        extends Modelica.Icons.Function;
      algorithm
        //  u := cv_const*(state.T - T0) - reference_p/d_const;
        u := cv_const*(state.T - T0);
        annotation (Documentation(info="<html>
<p>
This function computes the specific internal energy of the fluid, but neglects the (small) influence of the pressure term p/d.
</p>
</html>"));
      end specificInternalEnergy;

      redeclare function extends specificEntropy "Return specific entropy"
        extends Modelica.Icons.Function;
      algorithm
        s := cv_const*Modelica.Math.log(state.T/T0);
      end specificEntropy;

      redeclare function extends specificGibbsEnergy
        "Return specific Gibbs energy"
        extends Modelica.Icons.Function;
      algorithm
        g := specificEnthalpy(state) - state.T*specificEntropy(state);
      end specificGibbsEnergy;

      redeclare function extends specificHelmholtzEnergy
        "Return specific Helmholtz energy"
        extends Modelica.Icons.Function;
      algorithm
        f := specificInternalEnergy(state) - state.T*specificEntropy(state);
      end specificHelmholtzEnergy;

      redeclare function extends isentropicEnthalpy "Return isentropic enthalpy"
      algorithm
        h_is := cp_const*(temperature(refState) - T0);
      end isentropicEnthalpy;

      redeclare function extends isobaricExpansionCoefficient
        "Returns overall the isobaric expansion coefficient beta"
      algorithm
        beta := 0.0;
      end isobaricExpansionCoefficient;

      redeclare function extends isothermalCompressibility
        "Returns overall the isothermal compressibility factor"
      algorithm
        kappa := 0;
      end isothermalCompressibility;

      redeclare function extends density_derp_T
        "Returns the partial derivative of density with respect to pressure at constant temperature"
      algorithm
        ddpT := 0;
      end density_derp_T;

      redeclare function extends density_derT_p
        "Returns the partial derivative of density with respect to temperature at constant pressure"
      algorithm
        ddTp := 0;
      end density_derT_p;

      redeclare function extends density_derX
        "Returns the partial derivative of density with respect to mass fractions at constant pressure and temperature"
      algorithm
        dddX := fill(0, nX);
      end density_derX;

      redeclare function extends molarMass "Return the molar mass of the medium"
      algorithm
        MM := MM_const;
      end molarMass;
    end PartialSimpleMedium;

    package Choices "Types, constants to define menu choices"
      extends Modelica.Icons.Package;

      type IndependentVariables = enumeration(
          T "Temperature",
          pT "Pressure, Temperature",
          ph "Pressure, Specific Enthalpy",
          phX "Pressure, Specific Enthalpy, Mass Fraction",
          pTX "Pressure, Temperature, Mass Fractions",
          dTX "Density, Temperature, Mass Fractions")
        "Enumeration defining the independent variables of a medium";

      type ReferenceEnthalpy = enumeration(
          ZeroAt0K
            "The enthalpy is 0 at 0 K (default), if the enthalpy of formation is excluded",
          ZeroAt25C
            "The enthalpy is 0 at 25 degC, if the enthalpy of formation is excluded",
          UserDefined
            "The user-defined reference enthalpy is used at 293.15 K (25 degC)")
        "Enumeration defining the reference enthalpy of a medium" annotation (
          Evaluate=true);
      annotation (Documentation(info="<html>
<p>
Enumerations and data types for all types of fluids
</p>

<p>
Note: Reference enthalpy might have to be extended with enthalpy of formation.
</p>
</html>"));
    end Choices;

    package Types "Types to be used in fluid models"
      extends Modelica.Icons.Package;

      type AbsolutePressure = SI.AbsolutePressure (
          min=0,
          max=1.e8,
          nominal=1.e5,
          start=1.e5)
        "Type for absolute pressure with medium specific attributes";

      type Density = SI.Density (
          min=0,
          max=1.e5,
          nominal=1,
          start=1) "Type for density with medium specific attributes";

      type DynamicViscosity = SI.DynamicViscosity (
          min=0,
          max=1.e8,
          nominal=1.e-3,
          start=1.e-3)
        "Type for dynamic viscosity with medium specific attributes";

      type EnthalpyFlowRate = SI.EnthalpyFlowRate (
          nominal=1000.0,
          min=-1.0e8,
          max=1.e8) "Type for enthalpy flow rate with medium specific attributes";

      type MassFraction = Real (
          quantity="MassFraction",
          final unit="kg/kg",
          min=0,
          max=1,
          nominal=0.1) "Type for mass fraction with medium specific attributes";

      type MoleFraction = Real (
          quantity="MoleFraction",
          final unit="mol/mol",
          min=0,
          max=1,
          nominal=0.1) "Type for mole fraction with medium specific attributes";

      type MolarMass = SI.MolarMass (
          min=0.001,
          max=0.25,
          nominal=0.032) "Type for molar mass with medium specific attributes";

      type MolarVolume = SI.MolarVolume (
          min=1e-6,
          max=1.0e6,
          nominal=1.0) "Type for molar volume with medium specific attributes";

      type IsentropicExponent = SI.RatioOfSpecificHeatCapacities (
          min=1,
          max=500000,
          nominal=1.2,
          start=1.2)
        "Type for isentropic exponent with medium specific attributes";

      type SpecificEnergy = SI.SpecificEnergy (
          min=-1.0e8,
          max=1.e8,
          nominal=1.e6)
        "Type for specific energy with medium specific attributes";

      type SpecificInternalEnergy = SpecificEnergy
        "Type for specific internal energy with medium specific attributes";

      type SpecificEnthalpy = SI.SpecificEnthalpy (
          min=-1.0e10,
          max=1.e10,
          nominal=1.e6)
        "Type for specific enthalpy with medium specific attributes";

      type SpecificEntropy = SI.SpecificEntropy (
          min=-1.e7,
          max=1.e7,
          nominal=1.e3)
        "Type for specific entropy with medium specific attributes";

      type SpecificHeatCapacity = SI.SpecificHeatCapacity (
          min=0,
          max=1.e7,
          nominal=1.e3,
          start=1.e3)
        "Type for specific heat capacity with medium specific attributes";

      type Temperature = SI.Temperature (
          min=1,
          max=1.e4,
          nominal=300,
          start=288.15) "Type for temperature with medium specific attributes";

      type ThermalConductivity = SI.ThermalConductivity (
          min=0,
          max=500,
          nominal=1,
          start=1)
        "Type for thermal conductivity with medium specific attributes";

      type PrandtlNumber = SI.PrandtlNumber (
          min=1e-3,
          max=1e5,
          nominal=1.0) "Type for Prandtl number with medium specific attributes";

      type VelocityOfSound = SI.Velocity (
          min=0,
          max=1.e5,
          nominal=1000,
          start=1000)
        "Type for velocity of sound with medium specific attributes";

      type ExtraProperty = Real (min=0.0, start=1.0)
        "Type for unspecified, mass-specific property transported by flow";

      type ExtraPropertyFlowRate = Real (unit="kg/s")
        "Type for flow rate of unspecified, mass-specific property";

      type IsobaricExpansionCoefficient = Real (
          min=0,
          max=1.0e8,
          unit="1/K")
        "Type for isobaric expansion coefficient with medium specific attributes";

      type DipoleMoment = Real (
          min=0.0,
          max=2.0,
          unit="debye",
          quantity="ElectricDipoleMoment")
        "Type for dipole moment with medium specific attributes";

      type DerDensityByPressure = SI.DerDensityByPressure
        "Type for partial derivative of density with respect to pressure with medium specific attributes";

      type DerDensityByEnthalpy = SI.DerDensityByEnthalpy
        "Type for partial derivative of density with respect to enthalpy with medium specific attributes";

      type DerDensityByTemperature = SI.DerDensityByTemperature
        "Type for partial derivative of density with respect to temperature with medium specific attributes";

      package Basic "The most basic version of a record used in several degrees of detail"
        extends Icons.Package;

        record FluidConstants
          "Critical, triple, molecular and other standard data of fluid"
          extends Modelica.Icons.Record;
          String iupacName
            "Complete IUPAC name (or common name, if non-existent)";
          String casRegistryNumber
            "Chemical abstracts sequencing number (if it exists)";
          String chemicalFormula
            "Chemical formula, (brutto, nomenclature according to Hill";
          String structureFormula "Chemical structure formula";
          MolarMass molarMass "Molar mass";
        end FluidConstants;
      end Basic;

      package IdealGas "The ideal gas version of a record used in several degrees of detail"
        extends Icons.Package;

        record FluidConstants "Extended fluid constants"
          extends Modelica.Media.Interfaces.Types.Basic.FluidConstants;
          Temperature criticalTemperature "Critical temperature";
          AbsolutePressure criticalPressure "Critical pressure";
          MolarVolume criticalMolarVolume "Critical molar Volume";
          Real acentricFactor "Pitzer acentric factor";
          //   Temperature triplePointTemperature "Triple point temperature";
          //   AbsolutePressure triplePointPressure "Triple point pressure";
          Temperature meltingPoint "Melting point at 101325 Pa";
          Temperature normalBoilingPoint "Normal boiling point (at 101325 Pa)";
          DipoleMoment dipoleMoment
            "Dipole moment of molecule in Debye (1 debye = 3.33564e10-30 C.m)";
          Boolean hasIdealGasHeatCapacity=false
            "True if ideal gas heat capacity is available";
          Boolean hasCriticalData=false "True if critical data are known";
          Boolean hasDipoleMoment=false "True if a dipole moment known";
          Boolean hasFundamentalEquation=false "True if a fundamental equation";
          Boolean hasLiquidHeatCapacity=false
            "True if liquid heat capacity is available";
          Boolean hasSolidHeatCapacity=false
            "True if solid heat capacity is available";
          Boolean hasAccurateViscosityData=false
            "True if accurate data for a viscosity function is available";
          Boolean hasAccurateConductivityData=false
            "True if accurate data for thermal conductivity is available";
          Boolean hasVapourPressureCurve=false
            "True if vapour pressure data, e.g., Antoine coefficients are known";
          Boolean hasAcentricFactor=false
            "True if Pitzer acentric factor is known";
          SpecificEnthalpy HCRIT0=0.0
            "Critical specific enthalpy of the fundamental equation";
          SpecificEntropy SCRIT0=0.0
            "Critical specific entropy of the fundamental equation";
          SpecificEnthalpy deltah=0.0
            "Difference between specific enthalpy model (h_m) and f.eq. (h_f) (h_m - h_f)";
          SpecificEntropy deltas=0.0
            "Difference between specific enthalpy model (s_m) and f.eq. (s_f) (s_m - s_f)";
        end FluidConstants;
      end IdealGas;
    end Types;
    annotation (Documentation(info="<html>
<p>
This package provides basic interfaces definitions of media models for different
kind of media.
</p>
</html>"));
  end Interfaces;

  package Common "Data structures and fundamental functions for fluid properties"
    extends Modelica.Icons.Package;

    function smoothStep
      "Approximation of a general step, such that the characteristic is continuous and differentiable"
      extends Modelica.Icons.Function;
      input Real x "Abscissa value";
      input Real y1 "Ordinate value for x > 0";
      input Real y2 "Ordinate value for x < 0";
      input Real x_small(min=0) = 1e-5
        "Approximation of step for -x_small <= x <= x_small; x_small > 0 required";
      output Real y "Ordinate value to approximate y = if x > 0 then y1 else y2";
    algorithm
      y := smooth(1, if x > x_small then y1 else if x < -x_small then y2 else if
        abs(x_small) > 0 then (x/x_small)*((x/x_small)^2 - 3)*(y2 - y1)/4 + (y1
         + y2)/2 else (y1 + y2)/2);

      annotation (
        Inline=true,
        smoothOrder=1,
        Documentation(revisions="<html>
<ul>
<li><em>April 29, 2008</em>
    by <a href=\"mailto:Martin.Otter@DLR.de\">Martin Otter</a>:<br>
    Designed and implemented.</li>
<li><em>August 12, 2008</em>
    by <a href=\"mailto:Michael.Sielemann@dlr.de\">Michael Sielemann</a>:<br>
    Minor modification to cover the limit case <code>x_small -> 0</code> without division by zero.</li>
</ul>
</html>",   info="<html>
<p>
This function is used to approximate the equation
</p>
<blockquote><pre>
y = <strong>if</strong> x &gt; 0 <strong>then</strong> y1 <strong>else</strong> y2;
</pre></blockquote>

<p>
by a smooth characteristic, so that the expression is continuous and differentiable:
</p>

<blockquote><pre>
y = <strong>smooth</strong>(1, <strong>if</strong> x &gt;  x_small <strong>then</strong> y1 <strong>else</strong>
              <strong>if</strong> x &lt; -x_small <strong>then</strong> y2 <strong>else</strong> f(y1, y2));
</pre></blockquote>

<p>
In the region -x_small &lt; x &lt; x_small a 2nd order polynomial is used
for a smooth transition from y1 to y2.
</p>

<p>
If <strong>mass fractions</strong> X[:] are approximated with this function then this can be performed
for all <strong>nX</strong> mass fractions, instead of applying it for nX-1 mass fractions and computing
the last one by the mass fraction constraint sum(X)=1. The reason is that the approximating function has the
property that sum(X) = 1, provided sum(X_a) = sum(X_b) = 1
(and y1=X_a[i], y2=X_b[i]).
This can be shown by evaluating the approximating function in the abs(x) &lt; x_small
region (otherwise X is either X_a or X_b):
</p>

<blockquote><pre>
X[1]  = smoothStep(x, X_a[1] , X_b[1] , x_small);
X[2]  = smoothStep(x, X_a[2] , X_b[2] , x_small);
   ...
X[nX] = smoothStep(x, X_a[nX], X_b[nX], x_small);
</pre></blockquote>

<p>
or
</p>

<blockquote><pre>
X[1]  = c*(X_a[1]  - X_b[1])  + (X_a[1]  + X_b[1])/2
X[2]  = c*(X_a[2]  - X_b[2])  + (X_a[2]  + X_b[2])/2;
   ...
X[nX] = c*(X_a[nX] - X_b[nX]) + (X_a[nX] + X_b[nX])/2;
c     = (x/x_small)*((x/x_small)^2 - 3)/4
</pre></blockquote>

<p>
Summing all mass fractions together results in
</p>

<blockquote><pre>
sum(X) = c*(sum(X_a) - sum(X_b)) + (sum(X_a) + sum(X_b))/2
       = c*(1 - 1) + (1 + 1)/2
       = 1
</pre></blockquote>
</html>"));
    end smoothStep;
    annotation (Documentation(info="<html><h4>Package description</h4>
      <p>Package Modelica.Media.Common provides records and functions shared by many of the property sub-packages.
      High accuracy fluid property models share a lot of common structure, even if the actual models are different.
      Common data structures and computations shared by these property models are collected in this library.
   </p>

</html>",   revisions="<html>
      <ul>
      <li>First implemented: <em>July, 2000</em>
      by Hubertus Tummescheit
      for the ThermoFluid Library with help from Jonas Eborn and Falko Jens Wagner
      </li>
      <li>Code reorganization, enhanced documentation, additional functions: <em>December, 2002</em>
      by Hubertus Tummescheit and move to Modelica
                            properties library.</li>
      <li>Inclusion into Modelica.Media: September 2003</li>
      </ul>

      <address>Author: Hubertus Tummescheit,<br>
      Lund University<br>
      Department of Automatic Control<br>
      Box 118, 22100 Lund, Sweden<br>
      email: hubertus@control.lth.se
      </address>
</html>"));
  end Common;

    package Air "Medium models for air"
        extends Modelica.Icons.VariantsPackage;

      package MoistAir "Air: Moist air model (190 ... 647 K)"
        extends Interfaces.PartialCondensingGases(
          mediumName="Moist air",
          substanceNames={"water","air"},
          final reducedX=true,
          final singleState=false,
          reference_X={0.01,0.99},
          fluidConstants={IdealGases.Common.FluidData.H2O,IdealGases.Common.FluidData.N2},
          Temperature(min=190, max=647));

        import Modelica.Media.IdealGases.Common.Functions;
        constant Integer Water=1
          "Index of water (in substanceNames, massFractions X, etc.)";
        constant Integer Air=2
          "Index of air (in substanceNames, massFractions X, etc.)";
        //     constant SI.Pressure psat_low=saturationPressureWithoutLimits(200.0);
        //     constant SI.Pressure psat_high=saturationPressureWithoutLimits(422.16);
        constant Real k_mair=steam.MM/dryair.MM "Ratio of molar weights";

        constant IdealGases.Common.DataRecord dryair=IdealGases.Common.SingleGasesData.Air;
        constant IdealGases.Common.DataRecord steam=IdealGases.Common.SingleGasesData.H2O;
        constant SI.MolarMass[2] MMX={steam.MM,dryair.MM}
          "Molar masses of components";

        import Modelica.Media.Interfaces;
        import Modelica.Math;
        import Modelica.Constants;
        import Modelica.Media.IdealGases.Common.SingleGasNasa;
        import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;

        redeclare record extends ThermodynamicState
          "ThermodynamicState record for moist air"
        end ThermodynamicState;

        redeclare replaceable model extends BaseProperties(
          T(stateSelect=if preferredMediumStates then StateSelect.prefer else
                StateSelect.default),
          p(stateSelect=if preferredMediumStates then StateSelect.prefer else
                StateSelect.default),
          Xi(each stateSelect=if preferredMediumStates then StateSelect.prefer
                 else StateSelect.default),
          final standardOrderComponents=true) "Moist air base properties record"

          /* p, T, X = X[Water] are used as preferred states, since only then all
     other quantities can be computed in a recursive sequence.
     If other variables are selected as states, static state selection
     is no longer possible and non-linear algebraic equations occur.
      */
          MassFraction x_water "Mass of total water/mass of dry air";
          Real phi "Relative humidity";

      protected
          MassFraction X_liquid "Mass fraction of liquid or solid water";
          MassFraction X_steam "Mass fraction of steam water";
          MassFraction X_air "Mass fraction of air";
          MassFraction X_sat
            "Steam water mass fraction of saturation boundary in kg_water/kg_moistair";
          MassFraction x_sat
            "Steam water mass content of saturation boundary in kg_water/kg_dryair";
          AbsolutePressure p_steam_sat "Partial saturation pressure of steam";
        equation
          assert(T >= 190 and T <= 647, "
Temperature T is not in the allowed range
190.0 K <= (T ="     + String(T) + " K) <= 647.0 K
required from medium model \""     + mediumName + "\".");
          MM = 1/(Xi[Water]/MMX[Water] + (1.0 - Xi[Water])/MMX[Air]);

          p_steam_sat = min(saturationPressure(T), 0.999*p);
          X_sat = min(p_steam_sat*k_mair/max(100*Constants.eps, p - p_steam_sat)*(1
             - Xi[Water]), 1.0)
            "Water content at saturation with respect to actual water content";
          X_liquid = max(Xi[Water] - X_sat, 0.0);
          X_steam = Xi[Water] - X_liquid;
          X_air = 1 - Xi[Water];

          h = specificEnthalpy_pTX(
                p,
                T,
                Xi);
          R_s = dryair.R_s*(X_air/(1 - X_liquid)) + steam.R_s*X_steam/(1 - X_liquid);
          //
          u = h - R_s*T;
          d = p/(R_s*T);
          /* Note, u and d are computed under the assumption that the volume of the liquid
         water is negligible with respect to the volume of air and of steam
      */
          state.p = p;
          state.T = T;
          state.X = X;

          // these x are per unit mass of DRY air!
          x_sat = k_mair*p_steam_sat/max(100*Constants.eps, p - p_steam_sat);
          x_water = Xi[Water]/max(X_air, 100*Constants.eps);
          phi = p/p_steam_sat*Xi[Water]/(Xi[Water] + k_mair*X_air);
          annotation (Documentation(info="<html>
<p>This model computes thermodynamic properties of moist air from three independent (thermodynamic or/and numerical) state variables. Preferred numerical states are temperature T, pressure p and the reduced composition vector Xi, which contains the water mass fraction only. As an EOS the <strong>ideal gas law</strong> is used and associated restrictions apply. The model can also be used in the <strong>fog region</strong>, when moisture is present in its liquid state. However, it is assumed that the liquid water volume is negligible compared to that of the gas phase. Computation of thermal properties is based on property data of <a href=\"modelica://Modelica.Media.Air.DryAirNasa\"> dry air</a> and water (source: VDI-W&auml;rmeatlas), respectively. Besides the standard thermodynamic variables <strong>absolute and relative humidity</strong>, x_water and phi, respectively, are given by the model. Upper case X denotes absolute humidity with respect to mass of moist air while absolute humidity with respect to mass of dry air only is denoted by a lower case x throughout the model. See <a href=\"modelica://Modelica.Media.Air.MoistAir\">package description</a> for further information.</p>
</html>"));
        end BaseProperties;

        redeclare function setState_pTX
          "Return thermodynamic state as function of pressure p, temperature T and composition X"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input MassFraction X[:]=reference_X "Mass fractions";
          output ThermodynamicState state "Thermodynamic state";
        algorithm
          state := if size(X, 1) == nX then ThermodynamicState(
                p=p,
                T=T,
                X=X) else ThermodynamicState(
                p=p,
                T=T,
                X=cat(
                  1,
                  X,
                  {1 - sum(X)}));
          annotation (smoothOrder=2, Documentation(info="<html>
The <a href=\"modelica://Modelica.Media.Air.MoistAir.ThermodynamicState\">thermodynamic state record</a> is computed from pressure p, temperature T and composition X.
</html>"));
        end setState_pTX;

        redeclare function setState_phX
          "Return thermodynamic state as function of pressure p, specific enthalpy h and composition X"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input MassFraction X[:]=reference_X "Mass fractions";
          output ThermodynamicState state "Thermodynamic state";
        algorithm
          state := if size(X, 1) == nX then ThermodynamicState(
                p=p,
                T=T_phX(
                  p,
                  h,
                  X),
                X=X) else ThermodynamicState(
                p=p,
                T=T_phX(
                  p,
                  h,
                  X),
                X=cat(
                  1,
                  X,
                  {1 - sum(X)}));
          annotation (smoothOrder=2, Documentation(info="<html>
The <a href=\"modelica://Modelica.Media.Air.MoistAir.ThermodynamicState\">thermodynamic state record</a> is computed from pressure p, specific enthalpy h and composition X.
</html>"));
        end setState_phX;

        redeclare function setState_dTX
          "Return thermodynamic state as function of density d, temperature T and composition X"
          extends Modelica.Icons.Function;
          input Density d "Density";
          input Temperature T "Temperature";
          input MassFraction X[:]=reference_X "Mass fractions";
          output ThermodynamicState state "Thermodynamic state";
        algorithm
          state := if size(X, 1) == nX then ThermodynamicState(
                p=d*({steam.R_s,dryair.R_s}*X)*T,
                T=T,
                X=X) else ThermodynamicState(
                p=d*({steam.R_s,dryair.R_s}*cat(
                  1,
                  X,
                  {1 - sum(X)}))*T,
                T=T,
                X=cat(
                  1,
                  X,
                  {1 - sum(X)}));
          annotation (smoothOrder=2, Documentation(info="<html>
The <a href=\"modelica://Modelica.Media.Air.MoistAir.ThermodynamicState\">thermodynamic state record</a> is computed from density d, temperature T and composition X.
</html>"));
        end setState_dTX;

        redeclare function extends setSmoothState
          "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
        algorithm
          state := ThermodynamicState(
                p=Media.Common.smoothStep(
                  x,
                  state_a.p,
                  state_b.p,
                  x_small),
                T=Media.Common.smoothStep(
                  x,
                  state_a.T,
                  state_b.T,
                  x_small),
                X=Media.Common.smoothStep(
                  x,
                  state_a.X,
                  state_b.X,
                  x_small));
        end setSmoothState;

        function Xsaturation
          "Return absolute humidity per unit mass of moist air at saturation as a function of the thermodynamic state record"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output MassFraction X_sat "Steam mass fraction of sat. boundary";
        algorithm
          X_sat := k_mair/(state.p/min(saturationPressure(state.T), 0.999*state.p)
             - 1 + k_mair);
          annotation (smoothOrder=2, Documentation(info="<html>
Absolute humidity per unit mass of moist air at saturation is computed from pressure and temperature in the state record. Note, that unlike X_sat in the BaseProperties model this mass fraction refers to mass of moist air at saturation.
</html>"));
        end Xsaturation;

        function xsaturation
          "Return absolute humidity per unit mass of dry air at saturation as a function of the thermodynamic state record"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state record";
          output MassFraction x_sat "Absolute humidity per unit mass of dry air";
        algorithm
          x_sat := k_mair*saturationPressure(state.T)/max(100*Constants.eps, state.p
             - saturationPressure(state.T));
          annotation (smoothOrder=2, Documentation(info="<html>
Absolute humidity per unit mass of dry air at saturation is computed from pressure and temperature in the thermodynamic state record.
</html>"));
        end xsaturation;

        function xsaturation_pT
          "Return absolute humidity per unit mass of dry air at saturation as a function of pressure p and temperature T"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SI.Temperature T "Temperature";
          output MassFraction x_sat "Absolute humidity per unit mass of dry air";
        algorithm
          x_sat := k_mair*saturationPressure(T)/max(100*Constants.eps, p -
            saturationPressure(T));
          annotation (smoothOrder=2, Documentation(info="<html>
Absolute humidity per unit mass of dry air at saturation is computed from pressure and temperature.
</html>"));
        end xsaturation_pT;

        function massFraction_pTphi
          "Return steam mass fraction as a function of relative humidity phi and temperature T"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input Temperature T "Temperature";
          input Real phi "Relative humidity (0 ... 1.0)";
          output MassFraction X_steam "Absolute humidity, steam mass fraction";
      protected
          constant Real k=0.621964713077499 "Ratio of molar masses";
          AbsolutePressure psat=saturationPressure(T) "Saturation pressure";
        algorithm
          X_steam := phi*k/(k*phi + p/psat - phi);
          annotation (smoothOrder=2, Documentation(info="<html>
Absolute humidity per unit mass of moist air is computed from temperature, pressure and relative humidity.
</html>"));
        end massFraction_pTphi;

        function relativeHumidity_pTX
          "Return relative humidity as a function of pressure p, temperature T and composition X"
          extends Modelica.Icons.Function;
          input SI.Pressure p "Pressure";
          input SI.Temperature T "Temperature";
          input SI.MassFraction[:] X "Composition";
          output Real phi "Relative humidity";
      protected
          SI.Pressure p_steam_sat "Saturation pressure";
          SI.MassFraction X_air "Dry air mass fraction";
        algorithm
          p_steam_sat := min(saturationPressure(T), 0.999*p);
          X_air := 1 - X[Water];
          phi := max(0.0, min(1.0, p/p_steam_sat*X[Water]/(X[Water] + k_mair*X_air)));
          annotation (smoothOrder=2, Documentation(info="<html>
Relative humidity is computed from pressure, temperature and composition with 1.0 as the upper limit at saturation. Water mass fraction is the first entry in the composition vector.
</html>"));
        end relativeHumidity_pTX;

        function relativeHumidity
          "Return relative humidity as a function of the thermodynamic state record"
          extends Modelica.Icons.Function;
          input ThermodynamicState state "Thermodynamic state";
          output Real phi "Relative humidity";
        algorithm
          phi := relativeHumidity_pTX(
                state.p,
                state.T,
                state.X);
          annotation (smoothOrder=2, Documentation(info="<html>
Relative humidity is computed from the thermodynamic state record with 1.0 as the upper limit at saturation.
</html>"));
        end relativeHumidity;

        /*
    redeclare function setState_psX "Return thermodynamic state as function of p, s and composition X"
      extends Modelica.Icons.Function;
      input AbsolutePressure p "Pressure";
      input SpecificEntropy s "Specific entropy";
      input MassFraction X[:]=reference_X "Mass fractions";
      output ThermodynamicState state;
    algorithm
      state := if size(X,1) == nX then ThermodynamicState(p=p,T=T_psX(s,p,X),X=X)
        else ThermodynamicState(p=p,T=T_psX(p,s,X), X=cat(1,X,{1-sum(X)}));
    end setState_psX;
*/

        redeclare function extends gasConstant
          "Return ideal gas constant as a function from thermodynamic state, only valid for phi<1"

        algorithm
          R_s := dryair.R_s*(1 - state.X[Water]) + steam.R_s*state.X[Water];
          annotation (smoothOrder=2, Documentation(info="<html>
The ideal gas constant for moist air is computed from <a href=\"modelica://Modelica.Media.Air.MoistAir.ThermodynamicState\">thermodynamic state</a> assuming that all water is in the gas phase.
</html>"));
        end gasConstant;

        function gasConstant_X
          "Return ideal gas constant as a function from composition X"
          extends Modelica.Icons.Function;
          input SI.MassFraction X[:] "Gas phase composition";
          output SI.SpecificHeatCapacity R_s "Ideal gas constant";
        algorithm
          R_s := dryair.R_s*(1 - X[Water]) + steam.R_s*X[Water];
          annotation (smoothOrder=2, Documentation(info="<html>
The ideal gas constant for moist air is computed from the gas phase composition. The first entry in composition vector X is the steam mass fraction of the gas phase.
</html>"));
        end gasConstant_X;

        function saturationPressureLiquid
          "Return saturation pressure of water as a function of temperature T in the range of 273.16 to 647.096 K"

          extends Modelica.Icons.Function;
          input SI.Temperature Tsat "Saturation temperature";
          output SI.AbsolutePressure psat "Saturation pressure";
      protected
          SI.Temperature Tcritical=647.096 "Critical temperature";
          SI.AbsolutePressure pcritical=22.064e6 "Critical pressure";
          Real r1=(1 - Tsat/Tcritical) "Common subexpression";
          Real a[:]={-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,
              1.80122502} "Coefficients a[:]";
          Real n[:]={1.0,1.5,3.0,3.5,4.0,7.5} "Coefficients n[:]";
        algorithm
          psat := exp(((a[1]*r1^n[1] + a[2]*r1^n[2] + a[3]*r1^n[3] + a[4]*r1^n[4]
             + a[5]*r1^n[5] + a[6]*r1^n[6])*Tcritical)/Tsat)*pcritical;
          annotation (
            derivative=saturationPressureLiquid_der,
            Inline=false,
            smoothOrder=5,
            Documentation(info="<html>
<p>Saturation pressure of water above the triple point temperature is computed from temperature.</p>
<p>Source: A Saul, W Wagner: &quot;International equations for the saturation properties of ordinary water substance&quot;, equation 2.1</p>
</html>"));
        end saturationPressureLiquid;

        function saturationPressureLiquid_der
          "Derivative function for 'saturationPressureLiquid'"

          extends Modelica.Icons.Function;
          input SI.Temperature Tsat "Saturation temperature";
          input Real dTsat(unit="K/s") "Saturation temperature derivative";
          output Real psat_der(unit="Pa/s") "Saturation pressure derivative";
      protected
          SI.Temperature Tcritical=647.096 "Critical temperature";
          SI.AbsolutePressure pcritical=22.064e6 "Critical pressure";
          Real r1=(1 - Tsat/Tcritical) "Common subexpression 1";
          Real r1_der=-1/Tcritical*dTsat "Derivative of common subexpression 1";
          Real a[:]={-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,
              1.80122502} "Coefficients a[:]";
          Real n[:]={1.0,1.5,3.0,3.5,4.0,7.5} "Coefficients n[:]";
          Real r2=(a[1]*r1^n[1] + a[2]*r1^n[2] + a[3]*r1^n[3] + a[4]*r1^n[4] + a[5]
              *r1^n[5] + a[6]*r1^n[6]) "Common subexpression 2";
        algorithm
          // Approach used here is based on Baehr: "Thermodynamik", 12th edition p.204ff, "Method of Wagner"
          //psat := exp(((a[1]*r1^n[1] + a[2]*r1^n[2] + a[3]*r1^n[3] + a[4]*r1^n[4] + a[5]*r1^n[5] + a[6]*r1^n[6])*Tcritical)/Tsat) * pcritical;
          psat_der := exp((r2*Tcritical)/Tsat)*pcritical*((a[1]*(r1^(n[1] - 1)*n[1]
            *r1_der) + a[2]*(r1^(n[2] - 1)*n[2]*r1_der) + a[3]*(r1^(n[3] - 1)*n[3]*
            r1_der) + a[4]*(r1^(n[4] - 1)*n[4]*r1_der) + a[5]*(r1^(n[5] - 1)*n[5]*
            r1_der) + a[6]*(r1^(n[6] - 1)*n[6]*r1_der))*Tcritical/Tsat - r2*
            Tcritical*dTsat/Tsat^2);
          annotation (
            Inline=false,
            smoothOrder=5,
            Documentation(info="<html>
<p>Saturation pressure of water above the triple point temperature is computed from temperature.</p>
<p>Source: A Saul, W Wagner: &quot;International equations for the saturation properties of ordinary water substance&quot;, equation 2.1</p>
</html>"));
        end saturationPressureLiquid_der;

        function sublimationPressureIce
          "Return sublimation pressure of water as a function of temperature T between 190 and 273.16 K"

          extends Modelica.Icons.Function;
          input SI.Temperature Tsat "Sublimation temperature";
          output SI.AbsolutePressure psat "Sublimation pressure";
      protected
          SI.Temperature Ttriple=273.16 "Triple point temperature";
          SI.AbsolutePressure ptriple=611.657 "Triple point pressure";
          Real r1=Tsat/Ttriple "Common subexpression";
          Real a[:]={-13.9281690,34.7078238} "Coefficients a[:]";
          Real n[:]={-1.5,-1.25} "Coefficients n[:]";
        algorithm
          psat := exp(a[1] - a[1]*r1^n[1] + a[2] - a[2]*r1^n[2])*ptriple;
          annotation (
            Inline=false,
            smoothOrder=5,
            derivative=sublimationPressureIce_der,
            Documentation(info="<html>
<p>Sublimation pressure of water below the triple point temperature is computed from temperature.</p>
<p>Source: W Wagner, A Saul, A Pruss: &quot;International equations for the pressure along the melting and along the sublimation curve of ordinary water substance&quot;, equation 3.5</p>
</html>"));
        end sublimationPressureIce;

        function sublimationPressureIce_der
          "Derivative function for 'sublimationPressureIce'"

          extends Modelica.Icons.Function;
          input SI.Temperature Tsat "Sublimation temperature";
          input Real dTsat(unit="K/s") "Sublimation temperature derivative";
          output Real psat_der(unit="Pa/s") "Sublimation pressure derivative";
      protected
          SI.Temperature Ttriple=273.16 "Triple point temperature";
          SI.AbsolutePressure ptriple=611.657 "Triple point pressure";
          Real r1=Tsat/Ttriple "Common subexpression 1";
          Real r1_der=dTsat/Ttriple "Derivative of common subexpression 1";
          Real a[:]={-13.9281690,34.7078238} "Coefficients a[:]";
          Real n[:]={-1.5,-1.25} "Coefficients n[:]";
        algorithm
          //psat := exp(a[1] - a[1]*r1^n[1] + a[2] - a[2]*r1^n[2]) * ptriple;
          psat_der := exp(a[1] - a[1]*r1^n[1] + a[2] - a[2]*r1^n[2])*ptriple*(-(a[1]
            *(r1^(n[1] - 1)*n[1]*r1_der)) - (a[2]*(r1^(n[2] - 1)*n[2]*r1_der)));
          annotation (
            Inline=false,
            smoothOrder=5,
            Documentation(info="<html>
<p>Sublimation pressure of water below the triple point temperature is computed from temperature.</p>
<p>Source: W Wagner, A Saul, A Pruss: &quot;International equations for the pressure along the melting and along the sublimation curve of ordinary water substance&quot;, equation 3.5</p>
</html>"));
        end sublimationPressureIce_der;

        redeclare function extends saturationPressure
          "Return saturation pressure of water as a function of temperature T between 190 and 647.096 K"

        algorithm
          psat := Utilities.spliceFunction(
                saturationPressureLiquid(Tsat),
                sublimationPressureIce(Tsat),
                Tsat - 273.16,
                1.0);
          annotation (
            Inline=false,
            smoothOrder=5,
            derivative=saturationPressure_der,
            Documentation(info="<html>
Saturation pressure of water in the liquid and the solid region is computed using correlations. Functions for the
<a href=\"modelica://Modelica.Media.Air.MoistAir.sublimationPressureIce\">solid</a> and the <a href=\"modelica://Modelica.Media.Air.MoistAir.saturationPressureLiquid\"> liquid</a> region, respectively, are combined using the first derivative continuous <a href=\"modelica://Modelica.Media.Air.MoistAir.Utilities.spliceFunction\">spliceFunction</a>. This functions range of validity is from 190 to 647.096 K. For more information on the type of correlation used, see the documentation of the linked functions.
</html>"));
        end saturationPressure;

        function saturationPressure_der
          "Derivative function for 'saturationPressure'"
          extends Modelica.Icons.Function;
          input Temperature Tsat "Saturation temperature";
          input Real dTsat(unit="K/s") "Time derivative of saturation temperature";
          output Real psat_der(unit="Pa/s") "Saturation pressure";

        algorithm
          /*psat := Utilities.spliceFunction(saturationPressureLiquid(Tsat),sublimationPressureIce(Tsat),Tsat-273.16,1.0);*/
          psat_der := Utilities.spliceFunction_der(
                saturationPressureLiquid(Tsat),
                sublimationPressureIce(Tsat),
                Tsat - 273.16,
                1.0,
                saturationPressureLiquid_der(Tsat=Tsat, dTsat=dTsat),
                sublimationPressureIce_der(Tsat=Tsat, dTsat=dTsat),
                dTsat,
                0);
          annotation (
            Inline=false,
            smoothOrder=5,
            Documentation(info="<html>
Derivative function of <a href=\"modelica://Modelica.Media.Air.MoistAir.saturationPressure\">saturationPressure</a>
</html>"));
        end saturationPressure_der;

        function saturationTemperature
          "Return saturation temperature of water as a function of (partial) pressure p"
          extends Modelica.Icons.Function;
          input SI.Pressure p "Pressure";
          input SI.Temperature T_min=190 "Lower boundary of solution";
          input SI.Temperature T_max=647 "Upper boundary of solution";
          output SI.Temperature T "Saturation temperature";

      protected
          function f_nonlinear "Solve p(T) for T with given p"
            extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
            input SI.Pressure p "Pressure";
          algorithm
            y := saturationPressure(Tsat=u) - p;
          end f_nonlinear;

        algorithm
          T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
            function f_nonlinear(p=p), T_min, T_max);
          annotation (Documentation(info="<html>
Computes saturation temperature from (partial) pressure via numerical inversion of the function <a href=\"modelica://Modelica.Media.Air.MoistAir.saturationPressure\">saturationPressure</a>. Therefore additional inputs are required (or the defaults are used) for upper and lower temperature bounds.
</html>"));
        end saturationTemperature;

        redeclare function extends enthalpyOfVaporization
          "Return enthalpy of vaporization of water as a function of temperature T, 273.16 to 647.096 K"

      protected
          Real Tcritical=647.096 "Critical temperature";
          Real dcritical=322 "Critical density";
          Real pcritical=22.064e6 "Critical pressure";
          Real n[:]={1,1.5,3,3.5,4,7.5} "Powers in equation (1)";
          Real a[:]={-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,
              1.80122502} "Coefficients in equation (1) of [1]";
          Real m[:]={1/3,2/3,5/3,16/3,43/3,110/3} "Powers in equation (2)";
          Real b[:]={1.99274064,1.09965342,-0.510839303,-1.75493479,-45.5170352,-6.74694450e5}
            "Coefficients in equation (2) of [1]";
          Real o[:]={2/6,4/6,8/6,18/6,37/6,71/6} "Powers in equation (3)";
          Real c[:]={-2.03150240,-2.68302940,-5.38626492,-17.2991605,-44.7586581,-63.9201063}
            "Coefficients in equation (3) of [1]";
          Real tau=1 - T/Tcritical "Temperature expression";
          Real r1=(a[1]*Tcritical*tau^n[1])/T + (a[2]*Tcritical*tau^n[2])/T + (a[3]
              *Tcritical*tau^n[3])/T + (a[4]*Tcritical*tau^n[4])/T + (a[5]*
              Tcritical*tau^n[5])/T + (a[6]*Tcritical*tau^n[6])/T "Expression 1";
          Real r2=a[1]*n[1]*tau^n[1] + a[2]*n[2]*tau^n[2] + a[3]*n[3]*tau^n[3] + a[
              4]*n[4]*tau^n[4] + a[5]*n[5]*tau^n[5] + a[6]*n[6]*tau^n[6]
            "Expression 2";
          Real dp=dcritical*(1 + b[1]*tau^m[1] + b[2]*tau^m[2] + b[3]*tau^m[3] + b[
              4]*tau^m[4] + b[5]*tau^m[5] + b[6]*tau^m[6])
            "Density of saturated liquid";
          Real dpp=dcritical*exp(c[1]*tau^o[1] + c[2]*tau^o[2] + c[3]*tau^o[3] + c[
              4]*tau^o[4] + c[5]*tau^o[5] + c[6]*tau^o[6])
            "Density of saturated vapor";
        algorithm
          r0 := -(((dp - dpp)*exp(r1)*pcritical*(r2 + r1*tau))/(dp*dpp*tau))
            "Difference of equations (7) and (6)";
          annotation (
            smoothOrder=2,
            Documentation(info="<html>
<p>Enthalpy of vaporization of water is computed from temperature in the region of 273.16 to 647.096 K.</p>
<p>Source: W Wagner, A Pruss: \"International equations for the saturation properties of ordinary water substance. Revised according to the international temperature scale of 1990\" (1993).</p>
</html>"));
        end enthalpyOfVaporization;

        function HeatCapacityOfWater
          "Return specific heat capacity of water (liquid only) as a function of temperature T"
          extends Modelica.Icons.Function;
          input Temperature T "Temperature";
          output SpecificHeatCapacity cp_fl "Specific heat capacity of liquid";
        algorithm
          cp_fl := 1e3*(4.2166 - (T - 273.15)*(0.0033166 + (T - 273.15)*(0.00010295
             - (T - 273.15)*(1.3819e-6 + (T - 273.15)*7.3221e-9))));
          annotation (Documentation(info="<html>
The specific heat capacity of water (liquid and solid) is calculated using a
                 polynomial approach and data from VDI-Waermeatlas 8. Edition (Db1)
</html>"),     smoothOrder=2);
        end HeatCapacityOfWater;

        redeclare function extends enthalpyOfLiquid
          "Return enthalpy of liquid water as a function of temperature T(use enthalpyOfWater instead)"

        algorithm
          h := (T - 273.15)*1e3*(4.2166 - 0.5*(T - 273.15)*(0.0033166 + 0.333333*(T
             - 273.15)*(0.00010295 - 0.25*(T - 273.15)*(1.3819e-6 + 0.2*(T - 273.15)
            *7.3221e-9))));
          annotation (
            Inline=false,
            smoothOrder=5,
            Documentation(info="<html>
Specific enthalpy of liquid water is computed from temperature using a polynomial approach. Kept for compatibility reasons, better use <a href=\"modelica://Modelica.Media.Air.MoistAir.enthalpyOfWater\">enthalpyOfWater</a> instead.
</html>"));
        end enthalpyOfLiquid;

        redeclare function extends enthalpyOfGas
          "Return specific enthalpy of gas (air and steam) as a function of temperature T and composition X"

        algorithm
          h := Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=steam,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=46479.819 + 2501014.5)*X[Water] +
            Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=dryair,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=25104.684)*(1.0 - X[Water]);
          annotation (
            Inline=false,
            smoothOrder=5,
            Documentation(info="<html>
Specific enthalpy of moist air is computed from temperature, provided all water is in the gaseous state. The first entry in the composition vector X must be the mass fraction of steam. For a function that also covers the fog region please refer to <a href=\"modelica://Modelica.Media.Air.MoistAir.h_pTX\">h_pTX</a>.
</html>"));
        end enthalpyOfGas;

        redeclare function extends enthalpyOfCondensingGas
          "Return specific enthalpy of steam as a function of temperature T"

        algorithm
          h := Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=steam,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=46479.819 + 2501014.5);
          annotation (
            Inline=false,
            smoothOrder=5,
            Documentation(info="<html>
Specific enthalpy of steam is computed from temperature.
</html>"));
        end enthalpyOfCondensingGas;

        redeclare function extends enthalpyOfNonCondensingGas
          "Return specific enthalpy of dry air as a function of temperature T"

        algorithm
          h := Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=dryair,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=25104.684);
          annotation (
            Inline=false,
            smoothOrder=1,
            Documentation(info="<html>
Specific enthalpy of dry air is computed from temperature.
</html>"));
        end enthalpyOfNonCondensingGas;

        function enthalpyOfWater
          "Computes specific enthalpy of water (solid/liquid) near atmospheric pressure from temperature T"
          extends Modelica.Icons.Function;
          input SI.Temperature T "Temperature";
          output SI.SpecificEnthalpy h "Specific enthalpy of water";
        algorithm
          /*simple model assuming constant properties:
  heat capacity of liquid water:4200 J/kg
  heat capacity of solid water: 2050 J/kg
  enthalpy of fusion (liquid=>solid): 333000 J/kg*/

          h := Utilities.spliceFunction(
                4200*(T - 273.15),
                2050*(T - 273.15) - 333000,
                T - 273.16,
                0.1);
          annotation (derivative=enthalpyOfWater_der, Documentation(info="<html>
Specific enthalpy of water (liquid and solid) is computed from temperature using constant properties as follows:<br>
<ul>
<li>heat capacity of liquid water:4200 J/kg</li>
<li>heat capacity of solid water: 2050 J/kg</li>
<li>enthalpy of fusion (liquid=>solid): 333000 J/kg</li>
</ul>
Pressure is assumed to be around 1 bar. This function is usually used to determine the specific enthalpy of the liquid or solid fraction of moist air.
</html>"));
        end enthalpyOfWater;

        function enthalpyOfWater_der "Derivative function of enthalpyOfWater"
          extends Modelica.Icons.Function;
          input SI.Temperature T "Temperature";
          input Real dT(unit="K/s") "Time derivative of temperature";
          output Real dh(unit="J/(kg.s)") "Time derivative of specific enthalpy";
        algorithm
          /*simple model assuming constant properties:
  heat capacity of liquid water:4200 J/kg
  heat capacity of solid water: 2050 J/kg
  enthalpy of fusion (liquid=>solid): 333000 J/kg*/

          //h:=Utilities.spliceFunction(4200*(T-273.15),2050*(T-273.15)-333000,T-273.16,0.1);
          dh := Utilities.spliceFunction_der(
                4200*(T - 273.15),
                2050*(T - 273.15) - 333000,
                T - 273.16,
                0.1,
                4200*dT,
                2050*dT,
                dT,
                0);
          annotation (Documentation(info="<html>
Derivative function for <a href=\"modelica://Modelica.Media.Air.MoistAir.enthalpyOfWater\">enthalpyOfWater</a>.

</html>"));
        end enthalpyOfWater_der;

        redeclare function extends pressure
          "Returns pressure of ideal gas as a function of the thermodynamic state record"

        algorithm
          p := state.p;
          annotation (smoothOrder=2, Documentation(info="<html>
Pressure is returned from the thermodynamic state record input as a simple assignment.
</html>"));
        end pressure;

        redeclare function extends temperature
          "Return temperature of ideal gas as a function of the thermodynamic state record"

        algorithm
          T := state.T;
          annotation (smoothOrder=2, Documentation(info="<html>
Temperature is returned from the thermodynamic state record input as a simple assignment.
</html>"));
        end temperature;

        function T_phX
          "Return temperature as a function of pressure p, specific enthalpy h and composition X"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEnthalpy h "Specific enthalpy";
          input MassFraction[:] X "Mass fractions of composition";
          output Temperature T "Temperature";

      protected
          function f_nonlinear "Solve h_pTX(p,T,X) for T with given h"
            extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
            input AbsolutePressure p "Pressure";
            input SpecificEnthalpy h "Specific enthalpy";
            input MassFraction[:] X "Mass fractions of composition";
          algorithm
            y := h_pTX(p=p, T=u, X=X) - h;
          end f_nonlinear;

        algorithm
          T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
            function f_nonlinear(p=p, h=h, X=X[1:nXi]), 190, 647);
          annotation (Documentation(info="<html>
Temperature is computed from pressure, specific enthalpy and composition via numerical inversion of function <a href=\"modelica://Modelica.Media.Air.MoistAir.h_pTX\">h_pTX</a>.
</html>"));
        end T_phX;

        redeclare function extends density
          "Returns density of ideal gas as a function of the thermodynamic state record"

        algorithm
          d := state.p/(gasConstant(state)*state.T);
          annotation (smoothOrder=2, Documentation(info="<html>
Density is computed from pressure, temperature and composition in the thermodynamic state record applying the ideal gas law.
</html>"));
        end density;

        redeclare function extends specificEnthalpy
          "Return specific enthalpy of moist air as a function of the thermodynamic state record"

        algorithm
          h := h_pTX(
                state.p,
                state.T,
                state.X);
          annotation (smoothOrder=2, Documentation(info="<html>
Specific enthalpy of moist air is computed from the thermodynamic state record. The fog region is included for both, ice and liquid fog.
</html>"));
        end specificEnthalpy;

        function h_pTX
          "Return specific enthalpy of moist air as a function of pressure p, temperature T and composition X"
          extends Modelica.Icons.Function;
          input SI.Pressure p "Pressure";
          input SI.Temperature T "Temperature";
          input SI.MassFraction X[:] "Mass fractions of moist air";
          output SI.SpecificEnthalpy h "Specific enthalpy at p, T, X";
      protected
          SI.AbsolutePressure p_steam_sat "Partial saturation pressure of steam";
          SI.MassFraction X_sat "Absolute humidity per unit mass of moist air";
          SI.MassFraction X_liquid "Mass fraction of liquid water";
          SI.MassFraction X_steam "Mass fraction of steam water";
          SI.MassFraction X_air "Mass fraction of air";
        algorithm
          p_steam_sat := saturationPressure(T);
          //p_steam_sat :=min(saturationPressure(T), 0.999*p);
          X_sat := min(p_steam_sat*k_mair/max(100*Constants.eps, p - p_steam_sat)*(
            1 - X[Water]), 1.0);
          X_liquid := max(X[Water] - X_sat, 0.0);
          X_steam := X[Water] - X_liquid;
          X_air := 1 - X[Water];
          /* h        := {SingleGasNasa.h_Tlow(data=steam,  T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=46479.819+2501014.5),
               SingleGasNasa.h_Tlow(data=dryair, T=T, refChoice=ReferenceEnthalpy.UserDefined, h_off=25104.684)}*
    {X_steam, X_air} + enthalpyOfLiquid(T)*X_liquid;*/
          h := {Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=steam,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=46479.819 + 2501014.5),
            Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=dryair,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=25104.684)}*{X_steam,X_air} + enthalpyOfWater(T)*X_liquid;
          annotation (
            derivative=h_pTX_der,
            Inline=false,
            Documentation(info="<html>
Specific enthalpy of moist air is computed from pressure, temperature and composition with X[1] as the total water mass fraction. The fog region is included for both, ice and liquid fog.
</html>"));
        end h_pTX;

        function h_pTX_der "Derivative function of h_pTX"
          extends Modelica.Icons.Function;
          input SI.Pressure p "Pressure";
          input SI.Temperature T "Temperature";
          input SI.MassFraction X[:] "Mass fractions of moist air";
          input Real dp(unit="Pa/s") "Pressure derivative";
          input Real dT(unit="K/s") "Temperature derivative";
          input Real dX[:](each unit="1/s") "Composition derivative";
          output Real h_der(unit="J/(kg.s)") "Time derivative of specific enthalpy";
      protected
          SI.AbsolutePressure p_steam_sat "Partial saturation pressure of steam";
          SI.MassFraction X_sat "Absolute humidity per unit mass of moist air";
          SI.MassFraction X_liquid "Mass fraction of liquid water";
          SI.MassFraction X_steam "Mass fraction of steam water";
          SI.MassFraction X_air "Mass fraction of air";
          SI.MassFraction x_sat
            "Absolute humidity per unit mass of dry air at saturation";
          Real dX_steam(unit="1/s") "Time derivative of steam mass fraction";
          Real dX_air(unit="1/s") "Time derivative of dry air mass fraction";
          Real dX_liq(unit="1/s")
            "Time derivative of liquid/solid water mass fraction";
          Real dps(unit="Pa/s") "Time derivative of saturation pressure";
          Real dx_sat(unit="1/s")
            "Time derivative of absolute humidity per unit mass of dry air";
        algorithm
          p_steam_sat := saturationPressure(T);
          x_sat := p_steam_sat*k_mair/max(100*Modelica.Constants.eps, p -
            p_steam_sat);
          X_sat := min(x_sat*(1 - X[Water]), 1.0);
          X_liquid := Utilities.smoothMax(
                X[Water] - X_sat,
                0.0,
                1e-5);
          X_steam := X[Water] - X_liquid;
          X_air := 1 - X[Water];

          dX_air := -dX[Water];
          dps := saturationPressure_der(Tsat=T, dTsat=dT);
          dx_sat := k_mair*(dps*(p - p_steam_sat) - p_steam_sat*(dp - dps))/(p -
            p_steam_sat)/(p - p_steam_sat);
          dX_liq := Utilities.smoothMax_der(
                X[Water] - X_sat,
                0.0,
                1e-5,
                (1 + x_sat)*dX[Water] - (1 - X[Water])*dx_sat,
                0,
                0);
          dX_steam := dX[Water] - dX_liq;

          h_der := X_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow_der(
                data=steam,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=46479.819 + 2501014.5,
                dT=dT) + dX_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=steam,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=46479.819 + 2501014.5) + X_air*
            Modelica.Media.IdealGases.Common.Functions.h_Tlow_der(
                data=dryair,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=25104.684,
                dT=dT) + dX_air*Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=dryair,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=25104.684) + X_liquid*enthalpyOfWater_der(T=T, dT=dT) +
            dX_liq*enthalpyOfWater(T);

          annotation (
            Inline=false,
            smoothOrder=1,
            Documentation(info="<html>
Derivative function for <a href=\"modelica://Modelica.Media.Air.MoistAir.h_pTX\">h_pTX</a>.
</html>"));
        end h_pTX_der;

        redeclare function extends isentropicExponent
          "Return isentropic exponent (only for gas fraction!)"
        algorithm
          gamma := specificHeatCapacityCp(state)/specificHeatCapacityCv(state);
        end isentropicExponent;

        function isentropicEnthalpyApproximation
          "Approximate calculation of h_is from upstream properties, downstream pressure, gas part only"
          extends Modelica.Icons.Function;
          input AbsolutePressure p2 "Downstream pressure";
          input ThermodynamicState state "Thermodynamic state at upstream location";
          output SpecificEnthalpy h_is "Isentropic enthalpy";
      protected
          SpecificEnthalpy h "Specific enthalpy at upstream location";
          IsentropicExponent gamma=isentropicExponent(state) "Isentropic exponent";
      protected
          MassFraction[nX] X "Complete X-vector";
        algorithm
          X := state.X;
          //  X := if reducedX then cat(1,state.X,{1-sum(state.X)}) else state.X;
          h := {Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=steam,
                T=state.T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=46479.819 + 2501014.5),
            Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=dryair,
                T=state.T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=25104.684)}*X;

          h_is := h + gamma/(gamma - 1.0)*(state.T*gasConstant(state))*((p2/state.p)
            ^((gamma - 1)/gamma) - 1.0);
        end isentropicEnthalpyApproximation;

        redeclare function extends specificInternalEnergy
          "Return specific internal energy of moist air as a function of the thermodynamic state record"
          extends Modelica.Icons.Function;
        algorithm
          u := specificInternalEnergy_pTX(
                state.p,
                state.T,
                state.X);

          annotation (smoothOrder=2, Documentation(info="<html>
Specific internal energy is determined from the thermodynamic state record, assuming that the liquid or solid water volume is negligible.
</html>"));
        end specificInternalEnergy;

        function specificInternalEnergy_pTX
          "Return specific internal energy of moist air as a function of pressure p, temperature T and composition X"
          extends Modelica.Icons.Function;
          input SI.Pressure p "Pressure";
          input SI.Temperature T "Temperature";
          input SI.MassFraction X[:] "Mass fractions of moist air";
          output SI.SpecificInternalEnergy u "Specific internal energy";
      protected
          SI.AbsolutePressure p_steam_sat "Partial saturation pressure of steam";
          SI.MassFraction X_liquid "Mass fraction of liquid water";
          SI.MassFraction X_steam "Mass fraction of steam water";
          SI.MassFraction X_air "Mass fraction of air";
          SI.MassFraction X_sat "Absolute humidity per unit mass of moist air";
          SI.SpecificHeatCapacity R_gas "Ideal gas constant";
        algorithm
          p_steam_sat := saturationPressure(T);
          X_sat := min(p_steam_sat*k_mair/max(100*Constants.eps, p - p_steam_sat)*(
            1 - X[Water]), 1.0);
          X_liquid := max(X[Water] - X_sat, 0.0);
          X_steam := X[Water] - X_liquid;
          X_air := 1 - X[Water];
          R_gas := dryair.R_s*X_air/(1 - X_liquid) + steam.R_s*X_steam/(1 - X_liquid);
          u := X_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=steam,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=46479.819 + 2501014.5) + X_air*
            Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=dryair,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=25104.684) + enthalpyOfWater(T)*X_liquid - R_gas*T;

          annotation (derivative=specificInternalEnergy_pTX_der, Documentation(info=
                 "<html>
Specific internal energy is determined from pressure p, temperature T and composition X, assuming that the liquid or solid water volume is negligible.
</html>"));
        end specificInternalEnergy_pTX;

        function specificInternalEnergy_pTX_der
          "Derivative function for specificInternalEnergy_pTX"
          extends Modelica.Icons.Function;
          input SI.Pressure p "Pressure";
          input SI.Temperature T "Temperature";
          input SI.MassFraction X[:] "Mass fractions of moist air";
          input Real dp(unit="Pa/s") "Pressure derivative";
          input Real dT(unit="K/s") "Temperature derivative";
          input Real dX[:](each unit="1/s") "Mass fraction derivatives";
          output Real u_der(unit="J/(kg.s)") "Specific internal energy derivative";
      protected
          SI.AbsolutePressure p_steam_sat "Partial saturation pressure of steam";
          SI.MassFraction X_liquid "Mass fraction of liquid water";
          SI.MassFraction X_steam "Mass fraction of steam water";
          SI.MassFraction X_air "Mass fraction of air";
          SI.MassFraction X_sat "Absolute humidity per unit mass of moist air";
          SI.SpecificHeatCapacity R_gas "Ideal gas constant";

          SI.MassFraction x_sat
            "Absolute humidity per unit mass of dry air at saturation";
          Real dX_steam(unit="1/s") "Time derivative of steam mass fraction";
          Real dX_air(unit="1/s") "Time derivative of dry air mass fraction";
          Real dX_liq(unit="1/s")
            "Time derivative of liquid/solid water mass fraction";
          Real dps(unit="Pa/s") "Time derivative of saturation pressure";
          Real dx_sat(unit="1/s")
            "Time derivative of absolute humidity per unit mass of dry air";
          Real dR_gas(unit="J/(kg.K.s)") "Time derivative of ideal gas constant";
        algorithm
          p_steam_sat := saturationPressure(T);
          x_sat := p_steam_sat*k_mair/max(100*Modelica.Constants.eps, p -
            p_steam_sat);
          X_sat := min(x_sat*(1 - X[Water]), 1.0);
          X_liquid := Utilities.spliceFunction(
                X[Water] - X_sat,
                0.0,
                X[Water] - X_sat,
                1e-6);
          X_steam := X[Water] - X_liquid;
          X_air := 1 - X[Water];
          R_gas := steam.R_s*X_steam/(1 - X_liquid) + dryair.R_s*X_air/(1 - X_liquid);

          dX_air := -dX[Water];
          dps := saturationPressure_der(Tsat=T, dTsat=dT);
          dx_sat := k_mair*(dps*(p - p_steam_sat) - p_steam_sat*(dp - dps))/(p -
            p_steam_sat)/(p - p_steam_sat);
          dX_liq := Utilities.spliceFunction_der(
                X[Water] - X_sat,
                0.0,
                X[Water] - X_sat,
                1e-6,
                (1 + x_sat)*dX[Water] - (1 - X[Water])*dx_sat,
                0.0,
                (1 + x_sat)*dX[Water] - (1 - X[Water])*dx_sat,
                0.0);
          dX_steam := dX[Water] - dX_liq;
          dR_gas := (steam.R_s*(dX_steam*(1 - X_liquid) + dX_liq*X_steam) + dryair.R_s*
            (dX_air*(1 - X_liquid) + dX_liq*X_air))/(1 - X_liquid)/(1 - X_liquid);

          u_der := X_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow_der(
                data=steam,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=46479.819 + 2501014.5,
                dT=dT) + dX_steam*Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=steam,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=46479.819 + 2501014.5) + X_air*
            Modelica.Media.IdealGases.Common.Functions.h_Tlow_der(
                data=dryair,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=25104.684,
                dT=dT) + dX_air*Modelica.Media.IdealGases.Common.Functions.h_Tlow(
                data=dryair,
                T=T,
                refChoice=ReferenceEnthalpy.UserDefined,
                h_off=25104.684) + X_liquid*enthalpyOfWater_der(T=T, dT=dT) +
            dX_liq*enthalpyOfWater(T) - dR_gas*T - R_gas*dT;
          annotation (Documentation(info="<html>
Derivative function for <a href=\"modelica://Modelica.Media.Air.MoistAir.specificInternalEnergy_pTX\">specificInternalEnergy_pTX</a>.
</html>"));
        end specificInternalEnergy_pTX_der;

        redeclare function extends specificEntropy
          "Return specific entropy from thermodynamic state record, only valid for phi<1"

        algorithm
          s := s_pTX(
                state.p,
                state.T,
                state.X);
          annotation (
            Inline=false,
            smoothOrder=2,
            Documentation(info="<html>
Specific entropy is calculated from the thermodynamic state record, assuming ideal gas behavior and including entropy of mixing. Liquid or solid water is not taken into account, the entire water content X[1] is assumed to be in the vapor state (relative humidity below 1.0).
</html>"));
        end specificEntropy;

        redeclare function extends specificGibbsEnergy
          "Return specific Gibbs energy as a function of the thermodynamic state record, only valid for phi<1"
          extends Modelica.Icons.Function;
        algorithm
          g := h_pTX(
                state.p,
                state.T,
                state.X) - state.T*specificEntropy(state);
          annotation (smoothOrder=2, Documentation(info="<html>
The Gibbs Energy is computed from the thermodynamic state record for moist air with a water content below saturation.
</html>"));
        end specificGibbsEnergy;

        redeclare function extends specificHelmholtzEnergy
          "Return specific Helmholtz energy as a function of the thermodynamic state record, only valid for phi<1"
          extends Modelica.Icons.Function;
        algorithm
          f := h_pTX(
                state.p,
                state.T,
                state.X) - gasConstant(state)*state.T - state.T*specificEntropy(
            state);
          annotation (smoothOrder=2, Documentation(info="<html>
The Specific Helmholtz Energy is computed from the thermodynamic state record for moist air with a water content below saturation.
</html>"));
        end specificHelmholtzEnergy;

        redeclare function extends specificHeatCapacityCp
          "Return specific heat capacity at constant pressure as a function of the thermodynamic state record"

      protected
          Real dT(unit="s/K") = 1.0;
        algorithm
          cp := h_pTX_der(
                state.p,
                state.T,
                state.X,
                0.0,
                1.0,
                zeros(size(state.X, 1)))*dT "Definition of cp: dh/dT @ constant p";
          //      cp:= SingleGasNasa.cp_Tlow(dryair, state.T)*(1-state.X[Water])
          //        + SingleGasNasa.cp_Tlow(steam, state.T)*state.X[Water];
          annotation (
            Inline=false,
            smoothOrder=2,
            Documentation(info="<html>
The specific heat capacity at constant pressure <strong>cp</strong> is computed from temperature and composition for a mixture of steam (X[1]) and dry air. All water is assumed to be in the vapor state.
</html>"));
        end specificHeatCapacityCp;

        redeclare function extends specificHeatCapacityCv
          "Return specific heat capacity at constant volume as a function of the thermodynamic state record"

        algorithm
          cv := Modelica.Media.IdealGases.Common.Functions.cp_Tlow(dryair, state.T)
            *(1 - state.X[Water]) +
            Modelica.Media.IdealGases.Common.Functions.cp_Tlow(steam, state.T)*
            state.X[Water] - gasConstant(state);
          annotation (
            Inline=false,
            smoothOrder=2,
            Documentation(info="<html>
The specific heat capacity at constant density <strong>cv</strong> is computed from temperature and composition for a mixture of steam (X[1]) and dry air. All water is assumed to be in the vapor state.
</html>"));
        end specificHeatCapacityCv;

      redeclare function extends dynamicViscosity
          "Return dynamic viscosity as a function of the thermodynamic state record, valid from 123.15 K to 1273.15 K"

        import Modelica.Math.Polynomials;
      algorithm
        eta := 1e-6*Polynomials.evaluateWithRange(
            {9.7391102886305869E-15,-3.1353724870333906E-11,4.3004876595642225E-08,
            -3.8228016291758240E-05,5.0427874367180762E-02,1.7239260139242528E+01},
            Cv.to_degC(123.15),
            Cv.to_degC(1273.15),
            Cv.to_degC(state.T));
        annotation (smoothOrder=2, Documentation(info="<html>
<p>Dynamic viscosity is computed from temperature using a simple polynomial for dry air. Range of validity is from 123.15 K to 1273.15 K. The influence of pressure and moisture is neglected.</p>
<p>Source: VDI Waermeatlas, 8th edition.</p>
</html>"));
      end dynamicViscosity;

      redeclare function extends thermalConductivity
          "Return thermal conductivity as a function of the thermodynamic state record, valid from 123.15 K to 1273.15 K"
        import Modelica.Math.Polynomials;
      algorithm
        lambda := 1e-3*Polynomials.evaluateWithRange(
            {6.5691470817717812E-15,-3.4025961923050509E-11,5.3279284846303157E-08,
            -4.5340839289219472E-05,7.6129675309037664E-02,2.4169481088097051E+01},
            Cv.to_degC(123.15),
            Cv.to_degC(1273.15),
            Cv.to_degC(state.T));

        annotation (smoothOrder=2, Documentation(info="<html>
<p>Thermal conductivity is computed from temperature using a simple polynomial for dry air. Range of validity is from 123.15 K to 1273.15 K. The influence of pressure and moisture is neglected.</p>
<p>Source: VDI Waermeatlas, 8th edition.</p>
</html>"));
      end thermalConductivity;

        redeclare function extends velocityOfSound
        algorithm
          a := sqrt(isentropicExponent(state)*gasConstant(state)*temperature(state));
          annotation (Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end velocityOfSound;

        redeclare function extends isobaricExpansionCoefficient

        algorithm
          beta := 1/temperature(state);
          annotation (Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end isobaricExpansionCoefficient;

        redeclare function extends isothermalCompressibility

        algorithm
          kappa := 1/pressure(state);
          annotation (Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end isothermalCompressibility;

        redeclare function extends density_derp_h

        algorithm
          ddph := 1/(gasConstant(state)*temperature(state));

          annotation (Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end density_derp_h;

        redeclare function extends density_derh_p

        algorithm
          ddhp := -density(state)/(specificHeatCapacityCp(state)*temperature(state));
          annotation (Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end density_derh_p;

        redeclare function extends density_derp_T

        algorithm
          ddpT := 1/(gasConstant(state)*temperature(state));

          annotation (Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end density_derp_T;

        redeclare function extends density_derT_p

        algorithm
          ddTp := -density(state)/temperature(state);
          annotation (Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end density_derT_p;

        redeclare function extends density_derX

        algorithm
          dddX[Water] := - pressure(state)*(steam.R_s - dryair.R_s)/(((steam.R_s - dryair.R_s)
            *state.X[Water] + dryair.R_s)^2*temperature(state));
          dddX[Air] := - pressure(state)*(dryair.R_s - steam.R_s)/((steam.R_s + (dryair.R_s - steam.R_s)*
            state.X[Air])^2*temperature(state));

          annotation (Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
<p>2019-05-14        Stefan Wischhusen: Corrected derivatives.</p>
</html>"));
        end density_derX;

        redeclare function extends molarMass
        algorithm
          MM := Modelica.Constants.R/Modelica.Media.Air.MoistAir.gasConstant(state);
          annotation (Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end molarMass;

        function T_psX
          "Return temperature as a function of pressure p, specific entropy s and composition X"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          input MassFraction[:] X "Mass fractions of composition";
          output Temperature T "Temperature";

      protected
          function f_nonlinear "Solve s_pTX(p,T,X) for T with given s"
            extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
            input AbsolutePressure p "Pressure";
            input SpecificEntropy s "Specific entropy";
            input MassFraction[:] X "Mass fractions of composition";
          algorithm
            y := s_pTX(p=p, T=u, X=X) - s;
          end f_nonlinear;

        algorithm
          T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
            function f_nonlinear(p=p, s=s, X=X[1:nX]), 190, 647);
          annotation (Documentation(info="<html>
Temperature is computed from pressure, specific entropy and composition via numerical inversion of function <a href=\"modelica://Modelica.Media.Air.MoistAir.s_pTX\">s_pTX</a>.
</html>",     revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end T_psX;

        redeclare function extends setState_psX
        algorithm
          state := if size(X, 1) == nX then ThermodynamicState(
                p=p,
                T=T_psX(
                  p,
                  s,
                  X),
                X=X) else ThermodynamicState(
                p=p,
                T=T_psX(
                  p,
                  s,
                  X),
                X=cat(
                  1,
                  X,
                  {1 - sum(X)}));
          annotation (smoothOrder=2, Documentation(info="<html>
The <a href=\"modelica://Modelica.Media.Air.MoistAir.ThermodynamicState\">thermodynamic state record</a> is computed from pressure p, specific enthalpy h and composition X.
</html>",     revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end setState_psX;

        function s_pTX
          "Return specific entropy of moist air as a function of pressure p, temperature T and composition X (only valid for phi<1)"
          extends Modelica.Icons.Function;
          input SI.Pressure p "Pressure";
          input SI.Temperature T "Temperature";
          input SI.MassFraction X[:] "Mass fractions of moist air";
          output SI.SpecificEntropy s "Specific entropy at p, T, X";
      protected
          MoleFraction[2] Y=massToMoleFractions(X, {steam.MM,dryair.MM})
            "Molar fraction";

        algorithm
          s:= Modelica.Media.IdealGases.Common.Functions.s0_Tlow(dryair, T)*(1 - X[Water])
          + Modelica.Media.IdealGases.Common.Functions.s0_Tlow(steam, T)*X[Water]
          - Modelica.Constants.R*(Utilities.smoothMax(X[Water]/MMX[Water],0.0,1e-9)*Modelica.Math.log(max(Y[Water], Modelica.Constants.eps)*p/reference_p)
          + Utilities.smoothMax((1 - X[Water])/MMX[Air],0.0,1e-9)*Modelica.Math.log(max(Y[Air], Modelica.Constants.eps)*p/reference_p));
          annotation (
            derivative=s_pTX_der,
            Inline=false,
            Documentation(info="<html>
Specific entropy of moist air is computed from pressure, temperature and composition with X[1] as the total water mass fraction.
</html>",     revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
<p>2019-05-14        Stefan Wischhusen: Corrected calculation.</p>
<p>2019-09-10        Stefan Wischhusen: Corrected pressure influence (p &lt; p_ref).</p>
</html>"),  Icon(graphics={Text(
                  extent={{-100,100},{100,-100}},
                  textColor={255,127,0},
                  textString="f")}));
        end s_pTX;

        function s_pTX_der
          "Return specific entropy of moist air as a function of pressure p, temperature T and composition X (only valid for phi<1)"
          extends Modelica.Icons.Function;
          input SI.Pressure p "Pressure";
          input SI.Temperature T "Temperature";
          input SI.MassFraction X[:] "Mass fractions of moist air";
          input Real dp(unit="Pa/s") "Derivative of pressure";
          input Real dT(unit="K/s") "Derivative of temperature";
          input Real dX[nX](each unit="1/s") "Derivative of mass fractions";
          output Real ds(unit="J/(kg.K.s)") "Specific entropy at p, T, X";
      protected
          MoleFraction[2] Y=massToMoleFractions(X, {steam.MM,dryair.MM})
            "Molar fraction";
          MolarMass MM "Molar mass";

        algorithm
          MM := MMX[Water]*MMX[Air]/(X[Water]*MMX[Air] + X[Air]*MMX[Water]);

          ds := IdealGases.Common.Functions.s0_Tlow_der(
            dryair,
            T,
            dT)*(1 - X[Water]) + IdealGases.Common.Functions.s0_Tlow_der(
            steam,
            T,
            dT)*X[Water] + Modelica.Media.IdealGases.Common.Functions.s0_Tlow(dryair, T)*dX[Air] + Modelica.Media.IdealGases.Common.Functions.s0_Tlow(steam, T)*dX[Water] - Modelica.Constants.R*(1/MMX[Water]*(Utilities.smoothMax_der(
            X[Water],
            0.0,
            1e-9,
            dX[Water],
            0.0,
            0.0)*(Modelica.Math.log(max(Y[Water], Modelica.Constants.eps)*p/reference_p) + MM/MMX[Air]) + dp/p*Utilities.smoothMax(
            X[Water],
            0.0,
            1e-9)) + 1/MMX[Air]*(Utilities.smoothMax_der(
            X[Air],
            0.0,
            1e-9,
            dX[Air],
            0.0,
            0.0)*(Modelica.Math.log(max(Y[Air], Modelica.Constants.eps)*p/reference_p) + MM/MMX[Water]) + dp/p*Utilities.smoothMax(
            X[Air],
            0.0,
            1e-9)));

          annotation (
            Inline=false,
            smoothOrder=1,
            Documentation(info="<html>
Specific entropy of moist air is computed from pressure, temperature and composition with X[1] as the total water mass fraction.
</html>",     revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
<p>2019-05-14        Stefan Wischhusen: Corrected calculation.</p>
<p>2019-09-10        Stefan Wischhusen: Corrected pressure influence (p &lt; p_ref).</p>
</html>"),  Icon(graphics={Text(
                  extent={{-100,100},{100,-100}},
                  textColor={255,127,0},
                  textString="f")}));
        end s_pTX_der;

        redeclare function extends isentropicEnthalpy
          "Isentropic enthalpy (only valid for phi<1)"
          extends Modelica.Icons.Function;
        algorithm
          h_is := Modelica.Media.Air.MoistAir.h_pTX(
                p_downstream,
                Modelica.Media.Air.MoistAir.T_psX(
                  p_downstream,
                  Modelica.Media.Air.MoistAir.specificEntropy(refState),
                  refState.X),
                refState.X);

          annotation (Icon(graphics={Text(
                  extent={{-100,100},{100,-100}},
                  textColor={255,127,0},
                  textString="f")}), Documentation(revisions="<html>
<p>2012-01-12        Stefan Wischhusen: Initial Release.</p>
</html>"));
        end isentropicEnthalpy;

        package Utilities "Utility functions"
          extends Modelica.Icons.UtilitiesPackage;
          function spliceFunction "Spline interpolation of two functions"
            extends Modelica.Icons.Function;
            input Real pos "Returned value for x-deltax >= 0";
            input Real neg "Returned value for x+deltax <= 0";
            input Real x "Function argument";
            input Real deltax=1 "Region around x with spline interpolation";
            output Real out;
        protected
            Real scaledX;
            Real scaledX1;
            Real y;
          algorithm
            scaledX1 := x/deltax;
            scaledX := scaledX1*Modelica.Math.asin(1);
            if scaledX1 <= -0.999999999 then
              y := 0;
            elseif scaledX1 >= 0.999999999 then
              y := 1;
            else
              y := (Modelica.Math.tanh(Modelica.Math.tan(scaledX)) + 1)/2;
            end if;
            out := pos*y + (1 - y)*neg;
            annotation (derivative=spliceFunction_der);
          end spliceFunction;

          function spliceFunction_der "Derivative of spliceFunction"
            extends Modelica.Icons.Function;
            input Real pos;
            input Real neg;
            input Real x;
            input Real deltax=1;
            input Real dpos;
            input Real dneg;
            input Real dx;
            input Real ddeltax=0;
            output Real out;
        protected
            Real scaledX;
            Real scaledX1;
            Real dscaledX1;
            Real y;
          algorithm
            scaledX1 := x/deltax;
            scaledX := scaledX1*Modelica.Math.asin(1);
            dscaledX1 := (dx - scaledX1*ddeltax)/deltax;
            if scaledX1 <= -0.99999999999 then
              y := 0;
            elseif scaledX1 >= 0.9999999999 then
              y := 1;
            else
              y := (Modelica.Math.tanh(Modelica.Math.tan(scaledX)) + 1)/2;
            end if;
            out := dpos*y + (1 - y)*dneg;
            if (abs(scaledX1) < 1) then
              out := out + (pos - neg)*dscaledX1*Modelica.Math.asin(1)/2/(
                Modelica.Math.cosh(Modelica.Math.tan(scaledX))*Modelica.Math.cos(
                scaledX))^2;
            end if;
          end spliceFunction_der;

          function smoothMax
            extends Modelica.Icons.Function;
            import Modelica.Math;

            input Real x1 "First argument of smooth max operator";
            input Real x2 "Second argument of smooth max operator";
            input Real dx
              "Approximate difference between x1 and x2, below which regularization starts";
            output Real y "Result of smooth max operator";
          algorithm
            y := max(x1, x2) + Math.log((exp((4/dx)*(x1 - max(x1, x2)))) + (exp((4/
              dx)*(x2 - max(x1, x2)))))/(4/dx);
            annotation (smoothOrder=2, Documentation(info="<html>
<p>An implementation of Kreisselmeier Steinhauser smooth maximum</p>
</html>"));
          end smoothMax;

          function smoothMax_der
            extends Modelica.Icons.Function;

            import Modelica.Math.exp;
            import Modelica.Math.log;

            input Real x1 "First argument of smooth max operator";
            input Real x2 "Second argument of smooth max operator";
            input Real dx
              "Approximate difference between x1 and x2, below which regularization starts";
            input Real dx1;
            input Real dx2;
            input Real ddx;
            output Real dy "Derivative of smooth max operator";
          algorithm
            dy := (if x1 > x2 then dx1 else dx2) + 0.25*(((4*(dx1 - (if x1 > x2
               then dx1 else dx2))/dx - 4*(x1 - max(x1, x2))*ddx/dx^2)*exp(4*(x1 -
              max(x1, x2))/dx) + (4*(dx2 - (if x1 > x2 then dx1 else dx2))/dx - 4*(
              x2 - max(x1, x2))*ddx/dx^2)*exp(4*(x2 - max(x1, x2))/dx))*dx/(exp(4*(
              x1 - max(x1, x2))/dx) + exp(4*(x2 - max(x1, x2))/dx)) + log(exp(4*(x1
               - max(x1, x2))/dx) + exp(4*(x2 - max(x1, x2))/dx))*ddx);

            annotation (Documentation(info="<html>
<p>An implementation of Kreisselmeier Steinhauser smooth maximum</p>
</html>"));
          end smoothMax_der;
        end Utilities;

        annotation (Documentation(info="<html>
<h4>Thermodynamic Model</h4>
<p>This package provides a full thermodynamic model of moist air including the fog region and temperatures below zero degC.
The governing assumptions in this model are:</p>
<ul>
<li>the perfect gas law applies</li>
<li>water volume other than that of steam is neglected</li></ul>
<p>All extensive properties are expressed in terms of the total mass in order to comply with other media in this library. However, for moist air it is rather common to express the absolute humidity in terms of mass of dry air only, which has advantages when working with charts. In addition, care must be taken, when working with mass fractions with respect to total mass, that all properties refer to the same water content when being used in mathematical operations (which is always the case if based on dry air only). Therefore two absolute humidities are computed in the <strong>BaseProperties</strong> model: <strong>X</strong> denotes the absolute humidity in terms of the total mass while <strong>x</strong> denotes the absolute humidity per unit mass of dry air. In addition, the relative humidity <strong>phi</strong> is also computed.</p>
<p>At the triple point temperature of water of 0.01 &deg;C or 273.16 K and a relative humidity greater than 1 fog may be present as liquid and as ice resulting in a specific enthalpy somewhere between those of the two isotherms for solid and liquid fog, respectively. For numerical reasons a coexisting mixture of 50% solid and 50% liquid fog is assumed in the fog region at the triple point in this model.</p>

<h4>Range of validity</h4>
<p>From the assumptions mentioned above it follows that the <strong>pressure</strong> should be in the region around <strong>atmospheric</strong> conditions or below (a few bars may still be fine though). Additionally a very high water content at low temperatures would yield incorrect densities, because the volume of the liquid or solid phase would not be negligible anymore. The model does not provide information on limits for water drop size in the fog region or transport information for the actual condensation or evaporation process in combination with surfaces. All excess water which is not in its vapour state is assumed to be still present in the air regarding its energy but not in terms of its spatial extent.<br><br>
The thermodynamic model may be used for <strong>temperatures</strong> ranging from <strong>190 ... 647 K</strong>. This holds for all functions unless otherwise stated in their description. However, although the model works at temperatures above the saturation temperature it is questionable to use the term \"relative humidity\" in this region. Please note, that although several functions compute pure water properties, they are designed to be used within the moist air medium model where properties are dominated by air and steam in their vapor states, and not for pure liquid water applications.</p>

<h4>Transport Properties</h4>
<p>Several additional functions that are not needed to describe the thermodynamic system, but are required to model transport processes, like heat and mass transfer, may be called. They usually neglect the moisture influence unless otherwise stated.</p>

<h4>Application</h4>
<p>The model's main area of application is all processes that involve moist air cooling under near atmospheric pressure with possible moisture condensation. This is the case in all domestic and industrial air conditioning applications. Another large domain of moist air applications covers all processes that deal with dehydration of bulk material using air as a transport medium. Engineering tasks involving moist air are often performed (or at least visualized) by using charts that contain all relevant thermodynamic data for a moist air system. These so called psychrometric charts can be generated from the medium properties in this package. The model <a href=\"modelica://Modelica.Media.Examples.PsychrometricData\">PsychrometricData</a> may be used for this purpose in order to obtain data for figures like those below (the plotting itself is not part of the model though).</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Media/Air/Mollier.png\"><br>
<img src=\"modelica://Modelica/Resources/Images/Media/Air/PsycroChart.png\">
</p>

<p>
<strong>Legend:</strong> blue - constant specific enthalpy, red - constant temperature, black - constant relative humidity</p>

</html>"));
      end MoistAir;
      annotation (Documentation(info="<html>
  <p>This package contains different medium models for air:</p>
<ul>
<li><strong>SimpleAir</strong><br>
    Simple dry air medium in a limited temperature range.</li>
<li><strong>DryAirNasa</strong><br>
    Dry air as an ideal gas from Media.IdealGases.MixtureGases.Air.</li>
<li><strong>MoistAir</strong><br>
    Moist air as an ideal gas mixture of steam and dry air with fog below and above the triple point temperature.</li>
</ul>
</html>"));
    end Air;

    package IdealGases "Data and models of ideal gases (single, fixed and dynamic mixtures) from NASA source"
      extends Modelica.Icons.VariantsPackage;

      package Common "Common packages and data for the ideal gas models"
        extends Modelica.Icons.Package;

      record DataRecord
        "Coefficient data record for properties of ideal gases based on NASA source"
        extends Modelica.Icons.Record;
        String name "Name of ideal gas";
        SI.MolarMass MM "Molar mass";
        SI.SpecificEnthalpy Hf "Enthalpy of formation at 298.15K";
        SI.SpecificEnthalpy H0 "H0(298.15K) - H0(0K)";
        SI.Temperature Tlimit "Temperature limit between low and high data sets";
        Real alow[7] "Low temperature coefficients a";
        Real blow[2] "Low temperature constants b";
        Real ahigh[7] "High temperature coefficients a";
        Real bhigh[2] "High temperature constants b";
        SI.SpecificHeatCapacity R_s "Gas constant";
        annotation (Documentation(info="<html>
<p>
This data record contains the coefficients for the
ideal gas equations according to:
</p>
<blockquote>
  <p>McBride B.J., Zehe M.J., and Gordon S. (2002): <strong>NASA Glenn Coefficients
  for Calculating Thermodynamic Properties of Individual Species</strong>. NASA
  report TP-2002-211556</p>
</blockquote>
<p>
The equations have the following structure:
</p>
<div><img src=\"modelica://Modelica/Resources/Images/Media/IdealGases/Common/singleEquations.png\"></div>
<p>
The polynomials for h(T) and s0(T) are derived via integration from the one for cp(T)  and contain the integration constants b1, b2 that define the reference specific enthalpy and entropy. For entropy differences the reference pressure p0 is arbitrary, but not for absolute entropies. It is chosen as 1 standard atmosphere (101325 Pa).
</p>
<p>
For most gases, the region of validity is from 200 K to 6000 K.
The equations are split into two regions that are separated
by Tlimit (usually 1000 K). In both regions the gas is described
by the data above. The two branches are continuous and in most
gases also differentiable at Tlimit.
</p>
</html>"));
      end DataRecord;

        package Functions "Basic Functions for ideal gases: cp, h, s, thermal conductivity, viscosity"
          extends Modelica.Icons.FunctionsPackage;

          constant Boolean excludeEnthalpyOfFormation=true
            "If true, enthalpy of formation Hf is not included in specific enthalpy h";

          constant Modelica.Media.Interfaces.Choices.ReferenceEnthalpy referenceChoice=Modelica.Media.Interfaces.Choices.ReferenceEnthalpy.ZeroAt0K
            "Choice of reference enthalpy";

          constant Modelica.Media.Interfaces.Types.SpecificEnthalpy h_offset=0.0
            "User defined offset for reference enthalpy, if referenceChoice = UserDefined";

          constant Integer methodForThermalConductivity(min=1,max=2)=1;

          function cp_T
            "Compute specific heat capacity at constant pressure from temperature and gas data"
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            output SI.SpecificHeatCapacity cp "Specific heat capacity at temperature T";
          algorithm
            cp := smooth(0,if T < data.Tlimit then data.R_s*(1/(T*T)*(data.alow[1] + T*(
              data.alow[2] + T*(1.*data.alow[3] + T*(data.alow[4] + T*(data.alow[5] + T
              *(data.alow[6] + data.alow[7]*T))))))) else data.R_s*(1/(T*T)*(data.ahigh[1]
               + T*(data.ahigh[2] + T*(1.*data.ahigh[3] + T*(data.ahigh[4] + T*(data.
              ahigh[5] + T*(data.ahigh[6] + data.ahigh[7]*T))))))));
            annotation (Inline=true,smoothOrder=2);
          end cp_T;

          function cp_Tlow
            "Compute specific heat capacity at constant pressure, low T region"
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            output SI.SpecificHeatCapacity cp "Specific heat capacity at temperature T";
          algorithm
            cp := data.R_s*(1/(T*T)*(data.alow[1] + T*(
              data.alow[2] + T*(1.*data.alow[3] + T*(data.alow[4] + T*(data.alow[5] + T
              *(data.alow[6] + data.alow[7]*T)))))));
            annotation (Inline=false, derivative(zeroDerivative=data) = cp_Tlow_der);
          end cp_Tlow;

          function cp_Tlow_der
            "Compute derivative of specific heat capacity at constant pressure, low T region"
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            input Real dT(unit="K/s") "Temperature derivative";
            output Real cp_der(unit="J/(kg.K.s)") "Derivative of specific heat capacity";
          algorithm
            cp_der := dT*data.R_s/(T*T*T)*(-2*data.alow[1] + T*(
              -data.alow[2] + T*T*(data.alow[4] + T*(2.*data.alow[5] + T
              *(3.*data.alow[6] + 4.*data.alow[7]*T)))));
            annotation(smoothOrder=2);
          end cp_Tlow_der;

          function h_T "Compute specific enthalpy from temperature and gas data; reference is decided by the
    refChoice input, or by the referenceChoice package constant by default"
            import Modelica.Media.Interfaces.Choices;
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            input Boolean exclEnthForm=excludeEnthalpyOfFormation
              "If true, enthalpy of formation Hf is not included in specific enthalpy h";
            input Modelica.Media.Interfaces.Choices.ReferenceEnthalpy
                                            refChoice=referenceChoice
              "Choice of reference enthalpy";
            input SI.SpecificEnthalpy h_off=h_offset
              "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
            output SI.SpecificEnthalpy h "Specific enthalpy at temperature T";

          algorithm
            h := smooth(0,(if T < data.Tlimit then data.R_s*((-data.alow[1] + T*(data.
              blow[1] + data.alow[2]*Math.log(T) + T*(1.*data.alow[3] + T*(0.5*data.
              alow[4] + T*(1/3*data.alow[5] + T*(0.25*data.alow[6] + 0.2*data.alow[7]*T))))))
              /T) else data.R_s*((-data.ahigh[1] + T*(data.bhigh[1] + data.ahigh[2]*
              Math.log(T) + T*(1.*data.ahigh[3] + T*(0.5*data.ahigh[4] + T*(1/3*data.
              ahigh[5] + T*(0.25*data.ahigh[6] + 0.2*data.ahigh[7]*T))))))/T)) + (if
              exclEnthForm then -data.Hf else 0.0) + (if (refChoice
               == Choices.ReferenceEnthalpy.ZeroAt0K) then data.H0 else 0.0) + (if
              refChoice == Choices.ReferenceEnthalpy.UserDefined then h_off else
                    0.0));
            annotation (Inline=false,smoothOrder=2);
          end h_T;

          function h_Tlow "Compute specific enthalpy, low T region; reference is decided by the
    refChoice input, or by the referenceChoice package constant by default"
            import Modelica.Media.Interfaces.Choices;
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            input Boolean exclEnthForm=excludeEnthalpyOfFormation
              "If true, enthalpy of formation Hf is not included in specific enthalpy h";
            input Modelica.Media.Interfaces.Choices.ReferenceEnthalpy
                                            refChoice=referenceChoice
              "Choice of reference enthalpy";
            input SI.SpecificEnthalpy h_off=h_offset
              "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
            output SI.SpecificEnthalpy h "Specific enthalpy at temperature T";

          algorithm
            h := data.R_s*((-data.alow[1] + T*(data.
              blow[1] + data.alow[2]*Math.log(T) + T*(1.*data.alow[3] + T*(0.5*data.
              alow[4] + T*(1/3*data.alow[5] + T*(0.25*data.alow[6] + 0.2*data.alow[7]*T))))))
              /T) + (if
              exclEnthForm then -data.Hf else 0.0) + (if (refChoice
               == Choices.ReferenceEnthalpy.ZeroAt0K) then data.H0 else 0.0) + (if
              refChoice == Choices.ReferenceEnthalpy.UserDefined then h_off else
                    0.0);
            annotation(Inline=false,smoothOrder=2);
          end h_Tlow;

          function h_Tlow_der "Compute derivative of specific enthalpy, low T region; reference is decided by the
    refChoice input, or by the referenceChoice package constant by default"
            import Modelica.Media.Interfaces.Choices;
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            input Boolean exclEnthForm=excludeEnthalpyOfFormation
              "If true, enthalpy of formation Hf is not included in specific enthalpy h";
            input Modelica.Media.Interfaces.Choices.ReferenceEnthalpy
                                            refChoice=referenceChoice
              "Choice of reference enthalpy";
            input SI.SpecificEnthalpy h_off=h_offset
              "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
            input Real dT(unit="K/s") "Temperature derivative";
            output Real h_der(unit="J/(kg.s)")
              "Derivative of specific enthalpy at temperature T";
          algorithm
            h_der := dT*Modelica.Media.IdealGases.Common.Functions.cp_Tlow(
                                data,T);
            annotation(Inline=true,smoothOrder=2);
          end h_Tlow_der;

          function s0_T "Compute specific entropy from temperature and gas data"
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            output SI.SpecificEntropy s "Specific entropy at temperature T";
          algorithm
            s := if T < data.Tlimit then data.R_s*(data.blow[2] - 0.5*data.alow[
              1]/(T*T) - data.alow[2]/T + data.alow[3]*Math.log(T) + T*(
              data.alow[4] + T*(0.5*data.alow[5] + T*(1/3*data.alow[6] + 0.25*data.alow[
              7]*T)))) else data.R_s*(data.bhigh[2] - 0.5*data.ahigh[1]/(T*T) - data.
              ahigh[2]/T + data.ahigh[3]*Math.log(T) + T*(data.ahigh[4]
               + T*(0.5*data.ahigh[5] + T*(1/3*data.ahigh[6] + 0.25*data.ahigh[7]*T))));
            annotation (Inline=true, smoothOrder=2);
          end s0_T;

          function s0_Tlow "Compute specific entropy, low T region"
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            output SI.SpecificEntropy s "Specific entropy at temperature T";
          algorithm
            s := data.R_s*(data.blow[2] - 0.5*data.alow[
              1]/(T*T) - data.alow[2]/T + data.alow[3]*Math.log(T) + T*(
              data.alow[4] + T*(0.5*data.alow[5] + T*(1/3*data.alow[6] + 0.25*data.alow[
              7]*T))));
            annotation (Inline=true, smoothOrder=2);
          end s0_Tlow;

          function s0_Tlow_der "Compute derivative of specific entropy, low T region"
            extends Modelica.Icons.Function;
            input IdealGases.Common.DataRecord data "Ideal gas data";
            input SI.Temperature T "Temperature";
            input Real T_der(unit="K/s") "Temperature derivative";
            output Real s_der(unit="J/(kg.K.s)") "Derivative of specific entropy at temperature T";
          algorithm
            s_der := data.R_s*T_der*(data.alow[1]/(T*T*T) + data.alow[2]/(T*T) + data.alow[3]/T + data.alow[4] + T*(data.alow[5] + T*(data.alow[6] + T*data.alow[7])));
            annotation (Inline=true);
          end s0_Tlow_der;

          function dynamicViscosityLowPressure
            "Dynamic viscosity of low pressure gases"
            extends Modelica.Icons.Function;
            input SI.Temperature T "Gas temperature";
            input SI.Temperature Tc "Critical temperature of gas";
            input SI.MolarMass M "Molar mass of gas";
            input SI.MolarVolume Vc "Critical molar volume of gas";
            input Real w "Acentric factor of gas";
            input Modelica.Media.Interfaces.Types.DipoleMoment mu
              "Dipole moment of gas molecule";
            input Real k =  0.0 "Special correction for highly polar substances";
            output SI.DynamicViscosity eta "Dynamic viscosity of gas";
        protected
            parameter Real Const1_SI=40.785*10^(-9.5)
              "Constant in formula for eta converted to SI units";
            parameter Real Const2_SI=131.3/1000.0
              "Constant in formula for mur converted to SI units";
            Real mur=Const2_SI*mu/sqrt(Vc*Tc)
              "Dimensionless dipole moment of gas molecule";
            Real Fc=1 - 0.2756*w + 0.059035*mur^4 + k
              "Factor to account for molecular shape and polarities of gas";
            Real Tstar "Dimensionless temperature defined by equation below";
            Real Ov "Viscosity collision integral for the gas";

          algorithm
            Tstar := 1.2593*T/Tc;
            Ov := 1.16145*Tstar^(-0.14874) + 0.52487*Modelica.Math.exp(-0.7732*Tstar) + 2.16178*Modelica.Math.exp(-2.43787
              *Tstar);
            eta := Const1_SI*Fc*sqrt(M*T)/(Vc^(2/3)*Ov);
            annotation (smoothOrder=2,
                        Documentation(info="<html>
<p>
The used formula are based on the method of Chung et al (1984, 1988) referred to in ref [1] chapter 9.
The formula 9-4.10 is the one being used. The Formula is given in non-SI units, the following conversion constants were used to
transform the formula to SI units:
</p>

<ul>
<li> <strong>Const1_SI:</strong> The factor 10^(-9.5) =10^(-2.5)*1e-7 where the
     factor 10^(-2.5) originates from the conversion of g/mol->kg/mol + cm^3/mol->m^3/mol
      and the factor 1e-7 is due to conversion from microPoise->Pa.s.</li>
<li>  <strong>Const2_SI:</strong> The factor 1/3.335641e-27 = 1e-3/3.335641e-30
      where the factor 3.335641e-30 comes from debye->C.m and
      1e-3 is due to conversion from cm^3/mol->m^3/mol</li>
</ul>

<h4>References</h4>
<p>
[1] Bruce E. Poling, John E. Prausnitz, John P. O'Connell, \"The Properties of Gases and Liquids\" 5th Ed. Mc Graw Hill.
</p>

<h4>Author</h4>
<p>T. Skoglund, Lund, Sweden, 2004-08-31</p>

</html>"));
          end dynamicViscosityLowPressure;

          function thermalConductivityEstimate
            "Thermal conductivity of polyatomic gases (Eucken and Modified Eucken correlation)"
            extends Modelica.Icons.Function;
            input Modelica.Media.Interfaces.Types.SpecificHeatCapacity Cp
              "Constant pressure heat capacity";
            input Modelica.Media.Interfaces.Types.DynamicViscosity eta
              "Dynamic viscosity";
            input Integer method(min=1,max=2)=1
              "1: Eucken Method, 2: Modified Eucken Method";
            input IdealGases.Common.DataRecord data "Ideal gas data";
            output Modelica.Media.Interfaces.Types.ThermalConductivity lambda
              "Thermal conductivity [W/(m.k)]";
          algorithm
            lambda := if method == 1 then eta*(Cp - data.R_s + (9/4)*data.R_s)
                                     else eta*(Cp - data.R_s)*(1.32 + 1.77/((Cp/data.R_s) - 1.0));
            annotation (smoothOrder=2,
                        Documentation(info="<html>
<p>
This function provides two similar methods for estimating the
thermal conductivity of polyatomic gases.
The Eucken method (input method == 1) gives good results for low temperatures,
but it tends to give an underestimated value of the thermal conductivity
(lambda) at higher temperatures.<br>
The Modified Eucken method (input method == 2) gives good results for
high-temperatures, but it tends to give an overestimated value of the
thermal conductivity (lambda) at low temperatures.
</p>
</html>"));
          end thermalConductivityEstimate;
        end Functions;

      partial package SingleGasNasa "Medium model of an ideal gas based on NASA source"
        extends Interfaces.PartialPureSubstance(
           ThermoStates=Modelica.Media.Interfaces.Choices.IndependentVariables.pT,
           redeclare final record FluidConstants =
              Modelica.Media.Interfaces.Types.IdealGas.FluidConstants,
           mediumName=data.name,
           substanceNames={data.name},
           singleState=false,
           Temperature(min=200, max=6000, start=500, nominal=500),
           SpecificEnthalpy(start=if Functions.referenceChoice==ReferenceEnthalpy.ZeroAt0K then data.H0 else
              if Functions.referenceChoice==ReferenceEnthalpy.UserDefined then Functions.h_offset else 0, nominal=1.0e5),
           Density(start=10, nominal=10),
           AbsolutePressure(start=10e5, nominal=10e5));

        redeclare record extends ThermodynamicState
          "Thermodynamic state variables for ideal gases"
          AbsolutePressure p "Absolute pressure of medium";
          Temperature T "Temperature of medium";
        end ThermodynamicState;
        import Modelica.Math;
        import Modelica.Media.Interfaces.Choices.ReferenceEnthalpy;

        constant IdealGases.Common.DataRecord data
          "Data record of ideal gas substance";

        constant FluidConstants[nS] fluidConstants "Constant data for the fluid";

        redeclare model extends BaseProperties(
         T(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default),
         p(stateSelect=if preferredMediumStates then StateSelect.prefer else StateSelect.default))
          "Base properties of ideal gas medium"
        equation
          assert(T >= 200 and T <= 6000, "
Temperature T (= "       + String(T) + " K) is not in the allowed range
200 K <= T <= 6000 K required from medium model \""       + mediumName + "\".
");       MM = data.MM;
          R_s = data.R_s;
          h = Modelica.Media.IdealGases.Common.Functions.h_T(
                  data, T,
                  Modelica.Media.IdealGases.Common.Functions.excludeEnthalpyOfFormation,
                  Modelica.Media.IdealGases.Common.Functions.referenceChoice,
                  Modelica.Media.IdealGases.Common.Functions.h_offset);
          u = h - R_s*T;

          // Has to be written in the form d=f(p,T) in order that static
          // state selection for p and T is possible
          d = p/(R_s*T);
          // connect state with BaseProperties
          state.T = T;
          state.p = p;
        end BaseProperties;

          redeclare function setState_pTX
          "Return thermodynamic state as function of p, T and composition X"
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input Temperature T "Temperature";
            input MassFraction X[:]=reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := ThermodynamicState(p=p,T=T);
            annotation(Inline=true,smoothOrder=2);
          end setState_pTX;

          redeclare function setState_phX
          "Return thermodynamic state as function of p, h and composition X"
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input SpecificEnthalpy h "Specific enthalpy";
            input MassFraction X[:]=reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := ThermodynamicState(p=p,T=T_h(h));
            annotation(Inline=true,smoothOrder=2);
          end setState_phX;

          redeclare function setState_psX
          "Return thermodynamic state as function of p, s and composition X"
            extends Modelica.Icons.Function;
            input AbsolutePressure p "Pressure";
            input SpecificEntropy s "Specific entropy";
            input MassFraction X[:]=reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := ThermodynamicState(p=p,T=T_ps(p,s));
            annotation(Inline=true,smoothOrder=2);
          end setState_psX;

          redeclare function setState_dTX
          "Return thermodynamic state as function of d, T and composition X"
            extends Modelica.Icons.Function;
            input Density d "Density";
            input Temperature T "Temperature";
            input MassFraction X[:]=reference_X "Mass fractions";
            output ThermodynamicState state;
          algorithm
            state := ThermodynamicState(p=d*data.R_s*T,T=T);
            annotation(Inline=true,smoothOrder=2);
          end setState_dTX;

            redeclare function extends setSmoothState "Return thermodynamic state so that it smoothly approximates: if x > 0 then state_a else state_b"
            algorithm
              state := ThermodynamicState(p=Media.Common.smoothStep(x, state_a.p, state_b.p, x_small),
                                          T=Media.Common.smoothStep(x, state_a.T, state_b.T, x_small));
              annotation(Inline=true,smoothOrder=2);
            end setSmoothState;

        redeclare function extends pressure "Return pressure of ideal gas"
        algorithm
          p := state.p;
          annotation(Inline=true,smoothOrder=2);
        end pressure;

        redeclare function extends temperature "Return temperature of ideal gas"
        algorithm
          T := state.T;
          annotation(Inline=true,smoothOrder=2);
        end temperature;

        redeclare function extends density "Return density of ideal gas"
        algorithm
          d := state.p/(data.R_s*state.T);
          annotation(Inline=true,smoothOrder=2);
        end density;

        redeclare function extends specificEnthalpy "Return specific enthalpy"
          extends Modelica.Icons.Function;
        algorithm
          h := Modelica.Media.IdealGases.Common.Functions.h_T(
                   data,state.T);
          annotation(Inline=true,smoothOrder=2);
        end specificEnthalpy;

        redeclare function extends specificInternalEnergy
          "Return specific internal energy"
          extends Modelica.Icons.Function;
        algorithm
          u := Modelica.Media.IdealGases.Common.Functions.h_T(
                   data,state.T) - data.R_s*state.T;
          annotation(Inline=true,smoothOrder=2);
        end specificInternalEnergy;

        redeclare function extends specificEntropy "Return specific entropy"
          extends Modelica.Icons.Function;
        algorithm
          s := Modelica.Media.IdealGases.Common.Functions.s0_T(
                    data, state.T) - data.R_s*Modelica.Math.log(state.p/reference_p);
          annotation(Inline=true,smoothOrder=2);
        end specificEntropy;

        redeclare function extends specificGibbsEnergy "Return specific Gibbs energy"
          extends Modelica.Icons.Function;
        algorithm
          g := Modelica.Media.IdealGases.Common.Functions.h_T(
                   data,state.T) - state.T*specificEntropy(state);
          annotation(Inline=true,smoothOrder=2);
        end specificGibbsEnergy;

        redeclare function extends specificHelmholtzEnergy
          "Return specific Helmholtz energy"
          extends Modelica.Icons.Function;
        algorithm
          f := Modelica.Media.IdealGases.Common.Functions.h_T(
                   data,state.T) - data.R_s*state.T - state.T*specificEntropy(state);
          annotation(Inline=true,smoothOrder=2);
        end specificHelmholtzEnergy;

        redeclare function extends specificHeatCapacityCp
          "Return specific heat capacity at constant pressure"
        algorithm
          cp := Modelica.Media.IdealGases.Common.Functions.cp_T(
                     data, state.T);
          annotation(Inline=true,smoothOrder=2);
        end specificHeatCapacityCp;

        redeclare function extends specificHeatCapacityCv
          "Compute specific heat capacity at constant volume from temperature and gas data"
        algorithm
          cv := Modelica.Media.IdealGases.Common.Functions.cp_T(
                     data, state.T) - data.R_s;
          annotation(Inline=true,smoothOrder=2);
        end specificHeatCapacityCv;

        redeclare function extends isentropicExponent "Return isentropic exponent"
        algorithm
          gamma := specificHeatCapacityCp(state)/specificHeatCapacityCv(state);
          annotation(Inline=true,smoothOrder=2);
        end isentropicExponent;

        redeclare function extends velocityOfSound "Return velocity of sound"
          extends Modelica.Icons.Function;
        algorithm
          a := sqrt(max(0,data.R_s*state.T*Modelica.Media.IdealGases.Common.Functions.cp_T(
                                              data, state.T)/specificHeatCapacityCv(state)));
          annotation(Inline=true,smoothOrder=2);
        end velocityOfSound;

        function isentropicEnthalpyApproximation
          "Approximate method of calculating h_is from upstream properties and downstream pressure"
          extends Modelica.Icons.Function;
          input SI.Pressure p2 "Downstream pressure";
          input ThermodynamicState state "Properties at upstream location";
          input Boolean exclEnthForm=Functions.excludeEnthalpyOfFormation
            "If true, enthalpy of formation Hf is not included in specific enthalpy h";
          input ReferenceEnthalpy refChoice=Functions.referenceChoice
            "Choice of reference enthalpy";
          input SpecificEnthalpy h_off=Functions.h_offset
            "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
          output SI.SpecificEnthalpy h_is "Isentropic enthalpy";
      protected
          IsentropicExponent gamma =  isentropicExponent(state) "Isentropic exponent";
        algorithm
          h_is := Modelica.Media.IdealGases.Common.Functions.h_T(
                      data,state.T,exclEnthForm,refChoice,h_off) +
            gamma/(gamma - 1.0)*state.p/density(state)*((p2/state.p)^((gamma - 1)/gamma) - 1.0);
          annotation(Inline=true,smoothOrder=2);
        end isentropicEnthalpyApproximation;

        redeclare function extends isentropicEnthalpy "Return isentropic enthalpy"
        input Boolean exclEnthForm=Functions.excludeEnthalpyOfFormation
            "If true, enthalpy of formation Hf is not included in specific enthalpy h";
        input ReferenceEnthalpy refChoice=Functions.referenceChoice
            "Choice of reference enthalpy";
        input SpecificEnthalpy h_off=Functions.h_offset
            "User defined offset for reference enthalpy, if referenceChoice = UserDefined";
        algorithm
          h_is := isentropicEnthalpyApproximation(p_downstream,refState,exclEnthForm,refChoice,h_off);
          annotation(Inline=true,smoothOrder=2);
        end isentropicEnthalpy;

        redeclare function extends isobaricExpansionCoefficient
          "Returns overall the isobaric expansion coefficient beta"
        algorithm
          beta := 1/state.T;
          annotation(Inline=true,smoothOrder=2);
        end isobaricExpansionCoefficient;

        redeclare function extends isothermalCompressibility
          "Returns overall the isothermal compressibility factor"
        algorithm
          kappa := 1.0/state.p;
          annotation(Inline=true,smoothOrder=2);
        end isothermalCompressibility;

        redeclare function extends density_derp_T
          "Returns the partial derivative of density with respect to pressure at constant temperature"
        algorithm
          ddpT := 1/(state.T*data.R_s);
          annotation(Inline=true,smoothOrder=2);
        end density_derp_T;

        redeclare function extends density_derT_p
          "Returns the partial derivative of density with respect to temperature at constant pressure"
        algorithm
          ddTp := -state.p/(state.T*state.T*data.R_s);
          annotation(Inline=true,smoothOrder=2);
        end density_derT_p;

        redeclare function extends density_derX
          "Returns the partial derivative of density with respect to mass fractions at constant pressure and temperature"
        algorithm
          dddX := fill(0,nX);
          annotation(Inline=true,smoothOrder=2);
        end density_derX;

        redeclare replaceable function extends dynamicViscosity "Dynamic viscosity"
        algorithm
          assert(fluidConstants[1].hasCriticalData,
          "Failed to compute dynamicViscosity: For the species \"" + mediumName + "\" no critical data is available.");
          assert(fluidConstants[1].hasDipoleMoment,
          "Failed to compute dynamicViscosity: For the species \"" + mediumName + "\" no critical data is available.");
          eta := Modelica.Media.IdealGases.Common.Functions.dynamicViscosityLowPressure(
                                             state.T,
                             fluidConstants[1].criticalTemperature,
                             fluidConstants[1].molarMass,
                             fluidConstants[1].criticalMolarVolume,
                             fluidConstants[1].acentricFactor,
                             fluidConstants[1].dipoleMoment);
          annotation (smoothOrder=2);
        end dynamicViscosity;

        redeclare replaceable function extends thermalConductivity
          "Thermal conductivity of gas"
        //  input IdealGases.Common.DataRecord data "Ideal gas data";
          input Integer method=Functions.methodForThermalConductivity
            "1: Eucken Method, 2: Modified Eucken Method";
        algorithm
          assert(fluidConstants[1].hasCriticalData,
          "Failed to compute thermalConductivity: For the species \"" + mediumName + "\" no critical data is available.");
          lambda := Modelica.Media.IdealGases.Common.Functions.thermalConductivityEstimate(
                                                specificHeatCapacityCp(state),
            dynamicViscosity(state), method=method,data=data);
          annotation (smoothOrder=2);
        end thermalConductivity;

        redeclare function extends molarMass "Return the molar mass of the medium"
        algorithm
          MM := data.MM;
          annotation(Inline=true,smoothOrder=2);
        end molarMass;

        function T_h "Compute temperature from specific enthalpy"
          extends Modelica.Icons.Function;
          input SpecificEnthalpy h "Specific enthalpy";
          output Temperature T "Temperature";

      protected
          function f_nonlinear "Solve h(data,T) for T with given h (use only indirectly via temperature_phX)"
            extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
            input DataRecord data "Ideal gas data";
            input SpecificEnthalpy h "Specific enthalpy";
          algorithm
            y := Functions.h_T(data=data, T=u) - h;
          end f_nonlinear;

        algorithm
          T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
            function f_nonlinear(data=data, h=h), 200, 6000);
        end T_h;

        function T_ps "Compute temperature from pressure and specific entropy"
          extends Modelica.Icons.Function;
          input AbsolutePressure p "Pressure";
          input SpecificEntropy s "Specific entropy";
          output Temperature T "Temperature";

      protected
          function f_nonlinear "Solve s(data,T) for T with given s (use only indirectly via temperature_psX)"
            extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
            input DataRecord data "Ideal gas data";
            input AbsolutePressure p "Pressure";
            input SpecificEntropy s "Specific entropy";
          algorithm
            y := Functions.s0_T(data=data, T=u) - data.R_s*Modelica.Math.log(p/reference_p) - s;
          end f_nonlinear;

        algorithm
          T := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
            function f_nonlinear(data=data, p=p, s=s), 200, 6000);
        end T_ps;
        annotation (
          Documentation(info="<html>
<p>
This model calculates medium properties
for an ideal gas of a single substance, or for an ideal
gas consisting of several substances where the
mass fractions are fixed. Independent variables
are temperature <strong>T</strong> and pressure <strong>p</strong>.
Only density is a function of T and p. All other quantities
are solely a function of T. The properties
are valid in the range:
</p>
<blockquote><pre>
200 K &le; T &le; 6000 K
</pre></blockquote>
<p>
The following quantities are always computed:
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>Variable</strong></td>
      <td><strong>Unit</strong></td>
      <td><strong>Description</strong></td></tr>
  <tr><td>h</td>
      <td>J/kg</td>
      <td>specific enthalpy h = h(T)</td></tr>
  <tr><td>u</td>
      <td>J/kg</td>
      <td>specific internal energy u = u(T)</td></tr>
  <tr><td>d</td>
      <td>kg/m^3</td>
      <td>density d = d(p,T)</td></tr>
</table>
<p>
For the other variables, see the functions in
Modelica.Media.IdealGases.Common.SingleGasNasa.
Note, dynamic viscosity and thermal conductivity are only provided
for gases that use a data record from Modelica.Media.IdealGases.FluidData.
Currently these are the following gases:
</p>
<blockquote><pre>
Ar
C2H2_vinylidene
C2H4
C2H5OH
C2H6
C3H6_propylene
C3H7OH
C3H8
C4H8_1_butene
C4H9OH
C4H10_n_butane
C5H10_1_pentene
C5H12_n_pentane
C6H6
C6H12_1_hexene
C6H14_n_heptane
C7H14_1_heptene
C8H10_ethylbenz
CH3OH
CH4
CL2
CO
CO2
F2
H2
H2O
He
N2
N2O
NH3
NO
O2
SO2
SO3
</pre></blockquote>
<p>
<strong>Sources for model and literature:</strong><br>
Original Data: Computer program for calculation of complex chemical
equilibrium compositions and applications. Part 1: Analysis
Document ID: 19950013764 N (95N20180) File Series: NASA Technical Reports
Report Number: NASA-RP-1311  E-8017  NAS 1.61:1311
Authors: Gordon, Sanford (NASA Lewis Research Center)
 Mcbride, Bonnie J. (NASA Lewis Research Center)
Published: Oct 01, 1994.
</p>
<p><strong>Known limits of validity:</strong><br>
The data is valid for
temperatures between 200K and 6000K.  A few of the data sets for
monatomic gases have a discontinuous 1st derivative at 1000K, but
this never caused problems so far.
</p>
<p>
This model has been copied from the ThermoFluid library
and adapted to the Modelica.Media package.
</p>
</html>"));
      end SingleGasNasa;

        package FluidData "Critical data, dipole moments and related data"
          extends Modelica.Icons.Package;
          import Modelica.Media.Interfaces.PartialMixtureMedium;
          import Modelica.Media.IdealGases.Common.SingleGasesData;

          constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants N2(
                               chemicalFormula =        "N2",
                               iupacName =              "unknown",
                               structureFormula =       "unknown",
                               casRegistryNumber =      "7727-37-9",
                               meltingPoint =            63.15,
                               normalBoilingPoint =      77.35,
                               criticalTemperature =    126.20,
                               criticalPressure =        33.98e5,
                               criticalMolarVolume =     90.10e-6,
                               acentricFactor =           0.037,
                               dipoleMoment =             0.0,
                               molarMass =              SingleGasesData.N2.MM,
                               hasDipoleMoment =       true,
                               hasIdealGasHeatCapacity=true,
                               hasCriticalData =       true,
                               hasAcentricFactor =     true);

          constant Modelica.Media.Interfaces.Types.IdealGas.FluidConstants H2O(
                               chemicalFormula =        "H2O",
                               iupacName =              "oxidane",
                               structureFormula =       "H2O",
                               casRegistryNumber =      "7732-18-5",
                               meltingPoint =           273.15,
                               normalBoilingPoint =     373.124,
                               criticalTemperature =    647.096,
                               criticalPressure =       220.64e5,
                               criticalMolarVolume =     55.95e-6,
                               acentricFactor =           0.344,
                               dipoleMoment =             1.8,
                               molarMass =              SingleGasesData.H2O.MM,
                               hasDipoleMoment =       true,
                               hasIdealGasHeatCapacity=true,
                               hasCriticalData =       true,
                               hasAcentricFactor =     true);
          annotation (Documentation(info="<html>
<p>
This package contains FluidConstants data records for the following 37 gases
(see also the description in
<a href=\"modelica://Modelica.Media.IdealGases\">Modelica.Media.IdealGases</a>):
</p>
<blockquote><pre>
Argon             Methane          Methanol       Carbon Monoxide  Carbon Dioxide
Acetylene         Ethylene         Ethanol        Ethane           Propylene
Propane           1-Propanol       1-Butene       N-Butane         1-Pentene
N-Pentane         Benzene          1-Hexene       N-Hexane         1-Heptane
N-Heptane         Ethylbenzene     N-Octane       Chlorine         Fluorine
Hydrogen          Steam            Helium         Ammonia          Nitric Oxide
Nitrogen Dioxide  Nitrogen         Nitrous        Oxide            Neon Oxygen
Sulfur Dioxide    Sulfur Trioxide
</pre></blockquote>

</html>"));
        end FluidData;

        package SingleGasesData "Ideal gas data based on the NASA Glenn coefficients"
          extends Modelica.Icons.Package;

          constant IdealGases.Common.DataRecord Air(
            name="Air",
            MM=0.0289651159,
            Hf=-4333.833858403446,
            H0=298609.6803431054,
            Tlimit=1000,
            alow={10099.5016,-196.827561,5.00915511,-0.00576101373,1.06685993e-005,-7.94029797e-009,
                2.18523191e-012},
            blow={-176.796731,-3.921504225},
            ahigh={241521.443,-1257.8746,5.14455867,-0.000213854179,7.06522784e-008,-1.07148349e-011,
                6.57780015e-016},
            bhigh={6462.26319,-8.147411905},
            R_s=287.0512249529787);

          constant IdealGases.Common.DataRecord H2O(
            name="H2O",
            MM=0.01801528,
            Hf=-13423382.81725291,
            H0=549760.6476280135,
            Tlimit=1000,
            alow={-39479.6083,575.573102,0.931782653,0.00722271286,-7.34255737e-006,
                4.95504349e-009,-1.336933246e-012},
            blow={-33039.7431,17.24205775},
            ahigh={1034972.096,-2412.698562,4.64611078,0.002291998307,-6.836830479999999e-007,
                9.426468930000001e-011,-4.82238053e-015},
            bhigh={-13842.86509,-7.97814851},
            R_s=461.5233290850878);

          constant IdealGases.Common.DataRecord N2(
            name="N2",
            MM=0.0280134,
            Hf=0,
            H0=309498.4543111511,
            Tlimit=1000,
            alow={22103.71497,-381.846182,6.08273836,-0.00853091441,1.384646189e-005,-9.62579362e-009,
                2.519705809e-012},
            blow={710.846086,-10.76003744},
            ahigh={587712.406,-2239.249073,6.06694922,-0.00061396855,1.491806679e-007,-1.923105485e-011,
                1.061954386e-015},
            bhigh={12832.10415,-15.86640027},
            R_s=296.8033869505308);
          annotation (Documentation(info="<html>
<p>This package contains ideal gas models for the 1241 ideal gases from</p>
<blockquote>
  <p>McBride B.J., Zehe M.J., and Gordon S. (2002): <strong>NASA Glenn Coefficients
  for Calculating Thermodynamic Properties of Individual Species</strong>. NASA
  report TP-2002-211556</p>
</blockquote>

<blockquote><pre>
Ag        BaOH+           C2H4O_ethylen_o DF      In2I4    Nb      ScO2
Ag+       Ba_OH_2         CH3CHO_ethanal  DOCl    In2I6    Nb+     Sc2O
Ag-       BaS             CH3COOH         DO2     In2O     Nb-     Sc2O2
Air       Ba2             OHCH2COOH       DO2-    K        NbCl5   Si
Al        Be              C2H5            D2      K+       NbO     Si+
Al+       Be+             C2H5Br          D2+     K-       NbOCl3  Si-
Al-       Be++            C2H6            D2-     KAlF4    NbO2    SiBr
AlBr      BeBr            CH3N2CH3        D2O     KBO2     Ne      SiBr2
AlBr2     BeBr2           C2H5OH          D2O2    KBr      Ne+     SiBr3
AlBr3     BeCl            CH3OCH3         D2S     KCN      Ni      SiBr4
AlC       BeCl2           CH3O2CH3        e-      KCl      Ni+     SiC
AlC2      BeF             CCN             F       KF       Ni-     SiC2
AlCl      BeF2            CNC             F+      KH       NiCl    SiCl
AlCl+     BeH             OCCN            F-      KI       NiCl2   SiCl2
AlCl2     BeH+            C2N2            FCN     Kli      NiO     SiCl3
AlCl3     BeH2            C2O             FCO     KNO2     NiS     SiCl4
AlF       BeI             C3              FO      KNO3     O       SiF
AlF+      BeI2            C3H3_1_propynl  FO2_FOO KNa      O+      SiFCl
AlFCl     BeN             C3H3_2_propynl  FO2_OFO KO       O-      SiF2
AlFCl2    BeO             C3H4_allene     F2      KOH      OD      SiF3
AlF2      BeOH            C3H4_propyne    F2O     K2       OD-     SiF4
AlF2-     BeOH+           C3H4_cyclo      F2O2    K2+      OH      SiH
AlF2Cl    Be_OH_2         C3H5_allyl      FS2F    K2Br2    OH+     SiH+
AlF3      BeS             C3H6_propylene  Fe      K2CO3    OH-     SiHBr3
AlF4-     Be2             C3H6_cyclo      Fe+     K2C2N2   O2      SiHCl
AlH       Be2Cl4          C3H6O_propylox  Fe_CO_5 K2Cl2    O2+     SiHCl3
AlHCl     Be2F4           C3H6O_acetone   FeCl    K2F2     O2-     SiHF
AlHCl2    Be2O            C3H6O_propanal  FeCl2   K2I2     O3      SiHF3
AlHF      Be2OF2          C3H7_n_propyl   FeCl3   K2O      P       SiHI3
AlHFCl    Be2O2           C3H7_i_propyl   FeO     K2O+     P+      SiH2
AlHF2     Be3O3           C3H8            Fe_OH_2 K2O2     P-      SiH2Br2
AlH2      Be4O4           C3H8O_1propanol Fe2Cl4  K2O2H2   PCl     SiH2Cl2
AlH2Cl    Br              C3H8O_2propanol Fe2Cl6  K2SO4    PCl2    SiH2F2
AlH2F     Br+             CNCOCN          Ga      Kr       PCl2-   SiH2I2
AlH3      Br-             C3O2            Ga+     Kr+      PCl3    SiH3
AlI       BrCl            C4              GaBr    li       PCl5    SiH3Br
AlI2      BrF             C4H2_butadiyne  GaBr2   li+      PF      SiH3Cl
AlI3      BrF3            C4H4_1_3-cyclo  GaBr3   li-      PF+     SiH3F
AlN       BrF5            C4H6_butadiene  GaCl    liAlF4   PF-     SiH3I
AlO       BrO             C4H6_1butyne    GaCl2   liBO2    PFCl    SiH4
AlO+      OBrO            C4H6_2butyne    GaCl3   liBr     PFCl-   SiI
AlO-      BrOO            C4H6_cyclo      GaF     liCl     PFCl2   SiI2
AlOCl     BrO3            C4H8_1_butene   GaF2    liF      PFCl4   SiN
AlOCl2    Br2             C4H8_cis2_buten GaF3    liH      PF2     SiO
AlOF      BrBrO           C4H8_isobutene  GaH     liI      PF2-    SiO2
AlOF2     BrOBr           C4H8_cyclo      GaI     liN      PF2Cl   SiS
AlOF2-    C               C4H9_n_butyl    GaI2    liNO2    PF2Cl3  SiS2
AlOH      C+              C4H9_i_butyl    GaI3    liNO3    PF3     Si2
AlOHCl    C-              C4H9_s_butyl    GaO     liO      PF3Cl2  Si2C
AlOHCl2   CBr             C4H9_t_butyl    GaOH    liOF     PF4Cl   Si2F6
AlOHF     CBr2            C4H10_n_butane  Ga2Br2  liOH     PF5     Si2N
AlOHF2    CBr3            C4H10_isobutane Ga2Br4  liON     PH      Si3
AlO2      CBr4            C4N2            Ga2Br6  li2      PH2     Sn
AlO2-     CCl             C5              Ga2Cl2  li2+     PH2-    Sn+
Al_OH_2   CCl2            C5H6_1_3cyclo   Ga2Cl4  li2Br2   PH3     Sn-
Al_OH_2Cl CCl2Br2         C5H8_cyclo      Ga2Cl6  li2F2    PN      SnBr
Al_OH_2F  CCl3            C5H10_1_pentene Ga2F2   li2I2    PO      SnBr2
Al_OH_3   CCl3Br          C5H10_cyclo     Ga2F4   li2O     PO-     SnBr3
AlS       CCl4            C5H11_pentyl    Ga2F6   li2O+    POCl3   SnBr4
AlS2      CF              C5H11_t_pentyl  Ga2I2   li2O2    POFCl2  SnCl
Al2       CF+             C5H12_n_pentane Ga2I4   li2O2H2  POF2Cl  SnCl2
Al2Br6    CFBr3           C5H12_i_pentane Ga2I6   li2SO4   POF3    SnCl3
Al2C2     CFCl            CH3C_CH3_2CH3   Ga2O    li3+     PO2     SnCl4
Al2Cl6    CFClBr2         C6D5_phenyl     Ge      li3Br3   PO2-    SnF
Al2F6     CFCl2           C6D6            Ge+     li3Cl3   PS      SnF2
Al2I6     CFCl2Br         C6H2            Ge-     li3F3    P2      SnF3
Al2O      CFCl3           C6H5_phenyl     GeBr    li3I3    P2O3    SnF4
Al2O+     CF2             C6H5O_phenoxy   GeBr2   Mg       P2O4    SnI
Al2O2     CF2+            C6H6            GeBr3   Mg+      P2O5    SnI2
Al2O2+    CF2Br2          C6H5OH_phenol   GeBr4   MgBr     P3      SnI3
Al2O3     CF2Cl           C6H10_cyclo     GeCl    MgBr2    P3O6    SnI4
Al2S      CF2ClBr         C6H12_1_hexene  GeCl2   MgCl     P4      SnO
Al2S2     CF2Cl2          C6H12_cyclo     GeCl3   MgCl+    P4O6    SnO2
Ar        CF3             C6H13_n_hexyl   GeCl4   MgCl2    P4O7    SnS
Ar+       CF3+            C6H14_n_hexane  GeF     MgF      P4O8    SnS2
B         CF3Br           C7H7_benzyl     GeF2    MgF+     P4O9    Sn2
B+        CF3Cl           C7H8            GeF3    MgF2     P4O10   Sr
B-        CF4             C7H8O_cresol_mx GeF4    MgF2+    Pb      Sr+
BBr       CH+             C7H14_1_heptene GeH4    MgH      Pb+     SrBr
BBr2      CHBr3           C7H15_n_heptyl  GeI     MgI      Pb-     SrBr2
BBr3      CHCl            C7H16_n_heptane GeO     MgI2     PbBr    SrCl
BC        CHClBr2         C7H16_2_methylh GeO2    MgN      PbBr2   SrCl+
BC2       CHCl2           C8H8_styrene    GeS     MgO      PbBr3   SrCl2
BCl       CHCl2Br         C8H10_ethylbenz GeS2    MgOH     PbBr4   SrF
BCl+      CHCl3           C8H16_1_octene  Ge2     MgOH+    PbCl    SrF+
BClOH     CHF             C8H17_n_octyl   H       Mg_OH_2  PbCl2   SrF2
BCl_OH_2  CHFBr2          C8H18_n_octane  H+      MgS      PbCl3   SrH
BCl2      CHFCl           C8H18_isooctane H-      Mg2      PbCl4   SrI
BCl2+     CHFClBr         C9H19_n_nonyl   HAlO    Mg2F4    PbF     SrI2
BCl2OH    CHFCl2          C10H8_naphthale HAlO2   Mn       PbF2    SrO
BF        CHF2            C10H21_n_decyl  HBO     Mn+      PbF3    SrOH
BFCl      CHF2Br          C12H9_o_bipheny HBO+    Mo       PbF4    SrOH+
BFCl2     CHF2Cl          C12H10_biphenyl HBO2    Mo+      PbI     Sr_OH_2
BFOH      CHF3            Ca              HBS     Mo-      PbI2    SrS
BF_OH_2   CHI3            Ca+             HBS+    MoO      PbI3    Sr2
BF2       CH2             CaBr            HCN     MoO2     PbI4    Ta
BF2+      CH2Br2          CaBr2           HCO     MoO3     PbO     Ta+
BF2-      CH2Cl           CaCl            HCO+    MoO3-    PbO2    Ta-
BF2Cl     CH2ClBr         CaCl+           HCCN    Mo2O6    PbS     TaCl5
BF2OH     CH2Cl2          CaCl2           HCCO    Mo3O9    PbS2    TaO
BF3       CH2F            CaF             HCl     Mo4O12   Rb      TaO2
BF4-      CH2FBr          CaF+            HD      Mo5O15   Rb+     Ti
BH        CH2FCl          CaF2            HD+     N        Rb-     Ti+
BHCl      CH2F2           CaH             HDO     N+       RbBO2   Ti-
BHCl2     CH2I2           CaI             HDO2    N-       RbBr    TiCl
BHF       CH3             CaI2            HF      NCO      RbCl    TiCl2
BHFCl     CH3Br           CaO             HI      ND       RbF     TiCl3
BHF2      CH3Cl           CaO+            HNC     ND2      RbH     TiCl4
BH2       CH3F            CaOH            HNCO    ND3      RbI     TiO
BH2Cl     CH3I            CaOH+           HNO     NF       RbK     TiO+
BH2F      CH2OH           Ca_OH_2         HNO2    NF2      Rbli    TiOCl
BH3       CH2OH+          CaS             HNO3    NF3      RbNO2   TiOCl2
BH3NH3    CH3O            Ca2             HOCl    NH       RbNO3   TiO2
BH4       CH4             Cd              HOF     NH+      RbNa    U
BI        CH3OH           Cd+             HO2     NHF      RbO     UF
BI2       CH3OOH          Cl              HO2-    NHF2     RbOH    UF+
BI3       CI              Cl+             HPO     NH2      Rb2Br2  UF-
BN        CI2             Cl-             HSO3F   NH2F     Rb2Cl2  UF2
BO        CI3             ClCN            H2      NH3      Rb2F2   UF2+
BO-       CI4             ClF             H2+     NH2OH    Rb2I2   UF2-
BOCl      CN              ClF3            H2-     NH4+     Rb2O    UF3
BOCl2     CN+             ClF5            HBOH    NO       Rb2O2   UF3+
BOF       CN-             ClO             HCOOH   NOCl     Rb2O2H2 UF3-
BOF2      CNN             ClO2            H2F2    NOF      Rb2SO4  UF4
BOH       CO              Cl2             H2O     NOF3     Rn      UF4+
BO2       CO+             Cl2O            H2O+    NO2      Rn+     UF4-
BO2-      COCl            Co              H2O2    NO2-     S       UF5
B_OH_2    COCl2           Co+             H2S     NO2Cl    S+      UF5+
BS        COFCl           Co-             H2SO4   NO2F     S-      UF5-
BS2       COF2            Cr              H2BOH   NO3      SCl     UF6
B2        COHCl           Cr+             HB_OH_2 NO3-     SCl2    UF6-
B2C       COHF            Cr-             H3BO3   NO3F     SCl2+   UO
B2Cl4     COS             CrN             H3B3O3  N2       SD      UO+
B2F4      CO2             CrO             H3B3O6  N2+      SF      UOF
B2H       CO2+            CrO2            H3F3    N2-      SF+     UOF2
B2H2      COOH            CrO3            H3O+    NCN      SF-     UOF3
B2H3      CP              CrO3-           H4F4    N2D2_cis SF2     UOF4
B2H3_db   CS              Cs              H5F5    N2F2     SF2+    UO2
B2H4      CS2             Cs+             H6F6    N2F4     SF2-    UO2+
B2H4_db   C2              Cs-             H7F7    N2H2     SF3     UO2-
B2H5      C2+             CsBO2           He      NH2NO2   SF3+    UO2F
B2H5_db   C2-             CsBr            He+     N2H4     SF3-    UO2F2
B2H6      C2Cl            CsCl            Hg      N2O      SF4     UO3
B2O       C2Cl2           CsF             Hg+     N2O+     SF4+    UO3-
B2O2      C2Cl3           CsH             HgBr2   N2O3     SF4-    V
B2O3      C2Cl4           CsI             I       N2O4     SF5     V+
B2_OH_4   C2Cl6           Csli            I+      N2O5     SF5+    V-
B2S       C2F             CsNO2           I-      N3       SF5-    VCl4
B2S2      C2FCl           CsNO3           IF5     N3H      SF6     VN
B2S3      C2FCl3          CsNa            IF7     Na       SF6-    VO
B3H7_C2v  C2F2            CsO             I2      Na+      SH      VO2
B3H7_Cs   C2F2Cl2         CsOH            In      Na-      SH-     V4O10
B3H9      C2F3            CsRb            In+     NaAlF4   SN      W
B3N3H6    C2F3Cl          Cs2             InBr    NaBO2    SO      W+
B3O3Cl3   C2F4            Cs2Br2          InBr2   NaBr     SO-     W-
B3O3FCl2  C2F6            Cs2CO3          InBr3   NaCN     SOF2    WCl6
B3O3F2Cl  C2H             Cs2Cl2          InCl    NaCl     SO2     WO
B3O3F3    C2HCl           Cs2F2           InCl2   NaF      SO2-    WOCl4
B4H4      C2HCl3          Cs2I2           InCl3   NaH      SO2Cl2  WO2
B4H10     C2HF            Cs2O            InF     NaI      SO2FCl  WO2Cl2
B4H12     C2HFCl2         Cs2O+           InF2    Nali     SO2F2   WO3
B5H9      C2HF2Cl         Cs2O2           InF3    NaNO2    SO3     WO3-
Ba        C2HF3           Cs2O2H2         InH     NaNO3    S2      Xe
Ba+       C2H2_vinylidene Cs2SO4          InI     NaO      S2-     Xe+
BaBr      C2H2Cl2         Cu              InI2    NaOH     S2Cl2   Zn
BaBr2     C2H2FCl         Cu+             InI3    NaOH+    S2F2    Zn+
BaCl      C2H2F2          Cu-             InO     Na2      S2O     Zr
BaCl+     CH2CO_ketene    CuCl            InOH    Na2Br2   S3      Zr+
BaCl2     O_CH_2O         CuF             In2Br2  Na2Cl2   S4      Zr-
BaF       HO_CO_2OH       CuF2            In2Br4  Na2F2    S5      ZrN
BaF+      C2H3_vinyl      CuO             In2Br6  Na2I2    S6      ZrO
BaF2      CH2Br-COOH      Cu2             In2Cl2  Na2O     S7      ZrO+
BaH       C2H3Cl          Cu3Cl3          In2Cl4  Na2O+    S8      ZrO2
BaI       CH2Cl-COOH      D               In2Cl6  Na2O2    Sc
BaI2      C2H3F           D+              In2F2   Na2O2H2  Sc+
BaO       CH3CN           D-              In2F4   Na2SO4   Sc-
BaO+      CH3CO_acetyl    DBr             In2F6   Na3Cl3   ScO
BaOH      C2H4            DCl             In2I2   Na3F3    ScO+
</pre></blockquote>
</html>"));
        end SingleGasesData;
      annotation (Documentation(info="<html>
</html>"));
      end Common;
      annotation (Documentation(info="<html>
<p>This package contains data for the 1241 ideal gases from</p>
<blockquote>
  <p>McBride B.J., Zehe M.J., and Gordon S. (2002): <strong>NASA Glenn Coefficients
  for Calculating Thermodynamic Properties of Individual Species</strong>. NASA
  report TP-2002-211556</p>
</blockquote>
<p>Medium models for some of these gases are available in package
<a href=\"modelica://Modelica.Media.IdealGases.SingleGases\">IdealGases.SingleGases</a>
and some examples for mixtures are available in package <a href=\"modelica://Modelica.Media.IdealGases.MixtureGases\">IdealGases.MixtureGases</a>
</p>
<h4>Using and Adapting Medium Models</h4>
<p>
The data records allow computing the ideal gas specific enthalpy, specific entropy and heat capacity of the substances listed below. From them, even the Gibbs energy and equilibrium constants for reactions can be computed. Critical data that is needed for computing the viscosity and thermal conductivity is not included. In order to add mixtures or single substance medium packages that are
subtypes of
<a href=\"modelica://Modelica.Media.Interfaces.PartialMedium\">Interfaces.PartialMedium</a>
(i.e., can be utilized at all places where PartialMedium is defined),
a few additional steps have to be performed:
</p>
<ol>
<li>
All single gas media need to define a constant instance of record
<a href=\"modelica://Modelica.Media.Interfaces.PartialMedium.FluidConstants\">IdealGases.Common.SingleGasNasa.FluidConstants</a>.
For 37 ideal gases such records are provided in package
<a href=\"modelica://Modelica.Media.IdealGases.Common.FluidData\">IdealGases.Common.FluidData</a>.
For the other gases, such a record instance has to be provided by the user, e.g., by getting
the data from a commercial or public data base. A public source of the needed data is for example the <a href=\"http://webbook.nist.gov/chemistry/\"> NIST Chemistry WebBook</a></li>

<li>When the data is available, and a user has an instance of a
<a href=\"modelica://Modelica.Media.Interfaces.PartialMedium.FluidConstants\">FluidConstants</a> record filled with data, a medium package has to be written. Note that only the dipole moment, the acentric factor and critical data are necessary for the viscosity and thermal conductivity functions.</li>
<li><ul>
<li>For single components, a new package following the pattern in
<a href=\"modelica://Modelica.Media.IdealGases.SingleGases\">IdealGases.SingleGases</a> has to be created, pointing both to a data record for cp and to a user-defined fluidConstants record.</li>
<li>For mixtures of several components, a new package following the pattern in
<a href=\"modelica://Modelica.Media.IdealGases.MixtureGases\">IdealGases.MixtureGases</a> has to be created, building an array of data records for cp and an array of (partly) user-defined fluidConstants records.</li>
</ul></li>
</ol>
<p>Note that many properties can computed for the full set of 1241 gases listed below, but due to the missing viscosity and thermal conductivity functions, no fully Modelica.Media-compliant media can be defined.</p>
<p>
Data records for heat capacity, specific enthalpy and specific entropy exist for the following substances and ions:
</p>
<blockquote><pre>
Ag        BaOH+           C2H4O_ethylen_o DF      In2I4    Nb      ScO2
Ag+       Ba_OH_2         CH3CHO_ethanal  DOCl    In2I6    Nb+     Sc2O
Ag-       BaS             CH3COOH         DO2     In2O     Nb-     Sc2O2
Air       Ba2             OHCH2COOH       DO2-    K        NbCl5   Si
Al        Be              C2H5            D2      K+       NbO     Si+
Al+       Be+             C2H5Br          D2+     K-       NbOCl3  Si-
Al-       Be++            C2H6            D2-     KAlF4    NbO2    SiBr
AlBr      BeBr            CH3N2CH3        D2O     KBO2     Ne      SiBr2
AlBr2     BeBr2           C2H5OH          D2O2    KBr      Ne+     SiBr3
AlBr3     BeCl            CH3OCH3         D2S     KCN      Ni      SiBr4
AlC       BeCl2           CH3O2CH3        e-      KCl      Ni+     SiC
AlC2      BeF             CCN             F       KF       Ni-     SiC2
AlCl      BeF2            CNC             F+      KH       NiCl    SiCl
AlCl+     BeH             OCCN            F-      KI       NiCl2   SiCl2
AlCl2     BeH+            C2N2            FCN     Kli      NiO     SiCl3
AlCl3     BeH2            C2O             FCO     KNO2     NiS     SiCl4
AlF       BeI             C3              FO      KNO3     O       SiF
AlF+      BeI2            C3H3_1_propynl  FO2_FOO KNa      O+      SiFCl
AlFCl     BeN             C3H3_2_propynl  FO2_OFO KO       O-      SiF2
AlFCl2    BeO             C3H4_allene     F2      KOH      OD      SiF3
AlF2      BeOH            C3H4_propyne    F2O     K2       OD-     SiF4
AlF2-     BeOH+           C3H4_cyclo      F2O2    K2+      OH      SiH
AlF2Cl    Be_OH_2         C3H5_allyl      FS2F    K2Br2    OH+     SiH+
AlF3      BeS             C3H6_propylene  Fe      K2CO3    OH-     SiHBr3
AlF4-     Be2             C3H6_cyclo      Fe+     K2C2N2   O2      SiHCl
AlH       Be2Cl4          C3H6O_propylox  Fe_CO_5 K2Cl2    O2+     SiHCl3
AlHCl     Be2F4           C3H6O_acetone   FeCl    K2F2     O2-     SiHF
AlHCl2    Be2O            C3H6O_propanal  FeCl2   K2I2     O3      SiHF3
AlHF      Be2OF2          C3H7_n_propyl   FeCl3   K2O      P       SiHI3
AlHFCl    Be2O2           C3H7_i_propyl   FeO     K2O+     P+      SiH2
AlHF2     Be3O3           C3H8            Fe_OH_2 K2O2     P-      SiH2Br2
AlH2      Be4O4           C3H8O_1propanol Fe2Cl4  K2O2H2   PCl     SiH2Cl2
AlH2Cl    Br              C3H8O_2propanol Fe2Cl6  K2SO4    PCl2    SiH2F2
AlH2F     Br+             CNCOCN          Ga      Kr       PCl2-   SiH2I2
AlH3      Br-             C3O2            Ga+     Kr+      PCl3    SiH3
AlI       BrCl            C4              GaBr    li       PCl5    SiH3Br
AlI2      BrF             C4H2_butadiyne  GaBr2   li+      PF      SiH3Cl
AlI3      BrF3            C4H4_1_3-cyclo  GaBr3   li-      PF+     SiH3F
AlN       BrF5            C4H6_butadiene  GaCl    liAlF4   PF-     SiH3I
AlO       BrO             C4H6_1butyne    GaCl2   liBO2    PFCl    SiH4
AlO+      OBrO            C4H6_2butyne    GaCl3   liBr     PFCl-   SiI
AlO-      BrOO            C4H6_cyclo      GaF     liCl     PFCl2   SiI2
AlOCl     BrO3            C4H8_1_butene   GaF2    liF      PFCl4   SiN
AlOCl2    Br2             C4H8_cis2_buten GaF3    liH      PF2     SiO
AlOF      BrBrO           C4H8_isobutene  GaH     liI      PF2-    SiO2
AlOF2     BrOBr           C4H8_cyclo      GaI     liN      PF2Cl   SiS
AlOF2-    C               C4H9_n_butyl    GaI2    liNO2    PF2Cl3  SiS2
AlOH      C+              C4H9_i_butyl    GaI3    liNO3    PF3     Si2
AlOHCl    C-              C4H9_s_butyl    GaO     liO      PF3Cl2  Si2C
AlOHCl2   CBr             C4H9_t_butyl    GaOH    liOF     PF4Cl   Si2F6
AlOHF     CBr2            C4H10_n_butane  Ga2Br2  liOH     PF5     Si2N
AlOHF2    CBr3            C4H10_isobutane Ga2Br4  liON     PH      Si3
AlO2      CBr4            C4N2            Ga2Br6  li2      PH2     Sn
AlO2-     CCl             C5              Ga2Cl2  li2+     PH2-    Sn+
Al_OH_2   CCl2            C5H6_1_3cyclo   Ga2Cl4  li2Br2   PH3     Sn-
Al_OH_2Cl CCl2Br2         C5H8_cyclo      Ga2Cl6  li2F2    PN      SnBr
Al_OH_2F  CCl3            C5H10_1_pentene Ga2F2   li2I2    PO      SnBr2
Al_OH_3   CCl3Br          C5H10_cyclo     Ga2F4   li2O     PO-     SnBr3
AlS       CCl4            C5H11_pentyl    Ga2F6   li2O+    POCl3   SnBr4
AlS2      CF              C5H11_t_pentyl  Ga2I2   li2O2    POFCl2  SnCl
Al2       CF+             C5H12_n_pentane Ga2I4   li2O2H2  POF2Cl  SnCl2
Al2Br6    CFBr3           C5H12_i_pentane Ga2I6   li2SO4   POF3    SnCl3
Al2C2     CFCl            CH3C_CH3_2CH3   Ga2O    li3+     PO2     SnCl4
Al2Cl6    CFClBr2         C6D5_phenyl     Ge      li3Br3   PO2-    SnF
Al2F6     CFCl2           C6D6            Ge+     li3Cl3   PS      SnF2
Al2I6     CFCl2Br         C6H2            Ge-     li3F3    P2      SnF3
Al2O      CFCl3           C6H5_phenyl     GeBr    li3I3    P2O3    SnF4
Al2O+     CF2             C6H5O_phenoxy   GeBr2   Mg       P2O4    SnI
Al2O2     CF2+            C6H6            GeBr3   Mg+      P2O5    SnI2
Al2O2+    CF2Br2          C6H5OH_phenol   GeBr4   MgBr     P3      SnI3
Al2O3     CF2Cl           C6H10_cyclo     GeCl    MgBr2    P3O6    SnI4
Al2S      CF2ClBr         C6H12_1_hexene  GeCl2   MgCl     P4      SnO
Al2S2     CF2Cl2          C6H12_cyclo     GeCl3   MgCl+    P4O6    SnO2
Ar        CF3             C6H13_n_hexyl   GeCl4   MgCl2    P4O7    SnS
Ar+       CF3+            C6H14_n_hexane  GeF     MgF      P4O8    SnS2
B         CF3Br           C7H7_benzyl     GeF2    MgF+     P4O9    Sn2
B+        CF3Cl           C7H8            GeF3    MgF2     P4O10   Sr
B-        CF4             C7H8O_cresol_mx GeF4    MgF2+    Pb      Sr+
BBr       CH+             C7H14_1_heptene GeH4    MgH      Pb+     SrBr
BBr2      CHBr3           C7H15_n_heptyl  GeI     MgI      Pb-     SrBr2
BBr3      CHCl            C7H16_n_heptane GeO     MgI2     PbBr    SrCl
BC        CHClBr2         C7H16_2_methylh GeO2    MgN      PbBr2   SrCl+
BC2       CHCl2           C8H8_styrene    GeS     MgO      PbBr3   SrCl2
BCl       CHCl2Br         C8H10_ethylbenz GeS2    MgOH     PbBr4   SrF
BCl+      CHCl3           C8H16_1_octene  Ge2     MgOH+    PbCl    SrF+
BClOH     CHF             C8H17_n_octyl   H       Mg_OH_2  PbCl2   SrF2
BCl_OH_2  CHFBr2          C8H18_n_octane  H+      MgS      PbCl3   SrH
BCl2      CHFCl           C8H18_isooctane H-      Mg2      PbCl4   SrI
BCl2+     CHFClBr         C9H19_n_nonyl   HAlO    Mg2F4    PbF     SrI2
BCl2OH    CHFCl2          C10H8_naphthale HAlO2   Mn       PbF2    SrO
BF        CHF2            C10H21_n_decyl  HBO     Mn+      PbF3    SrOH
BFCl      CHF2Br          C12H9_o_bipheny HBO+    Mo       PbF4    SrOH+
BFCl2     CHF2Cl          C12H10_biphenyl HBO2    Mo+      PbI     Sr_OH_2
BFOH      CHF3            Ca              HBS     Mo-      PbI2    SrS
BF_OH_2   CHI3            Ca+             HBS+    MoO      PbI3    Sr2
BF2       CH2             CaBr            HCN     MoO2     PbI4    Ta
BF2+      CH2Br2          CaBr2           HCO     MoO3     PbO     Ta+
BF2-      CH2Cl           CaCl            HCO+    MoO3-    PbO2    Ta-
BF2Cl     CH2ClBr         CaCl+           HCCN    Mo2O6    PbS     TaCl5
BF2OH     CH2Cl2          CaCl2           HCCO    Mo3O9    PbS2    TaO
BF3       CH2F            CaF             HCl     Mo4O12   Rb      TaO2
BF4-      CH2FBr          CaF+            HD      Mo5O15   Rb+     Ti
BH        CH2FCl          CaF2            HD+     N        Rb-     Ti+
BHCl      CH2F2           CaH             HDO     N+       RbBO2   Ti-
BHCl2     CH2I2           CaI             HDO2    N-       RbBr    TiCl
BHF       CH3             CaI2            HF      NCO      RbCl    TiCl2
BHFCl     CH3Br           CaO             HI      ND       RbF     TiCl3
BHF2      CH3Cl           CaO+            HNC     ND2      RbH     TiCl4
BH2       CH3F            CaOH            HNCO    ND3      RbI     TiO
BH2Cl     CH3I            CaOH+           HNO     NF       RbK     TiO+
BH2F      CH2OH           Ca_OH_2         HNO2    NF2      Rbli    TiOCl
BH3       CH2OH+          CaS             HNO3    NF3      RbNO2   TiOCl2
BH3NH3    CH3O            Ca2             HOCl    NH       RbNO3   TiO2
BH4       CH4             Cd              HOF     NH+      RbNa    U
BI        CH3OH           Cd+             HO2     NHF      RbO     UF
BI2       CH3OOH          Cl              HO2-    NHF2     RbOH    UF+
BI3       CI              Cl+             HPO     NH2      Rb2Br2  UF-
BN        CI2             Cl-             HSO3F   NH2F     Rb2Cl2  UF2
BO        CI3             ClCN            H2      NH3      Rb2F2   UF2+
BO-       CI4             ClF             H2+     NH2OH    Rb2I2   UF2-
BOCl      CN              ClF3            H2-     NH4+     Rb2O    UF3
BOCl2     CN+             ClF5            HBOH    NO       Rb2O2   UF3+
BOF       CN-             ClO             HCOOH   NOCl     Rb2O2H2 UF3-
BOF2      CNN             ClO2            H2F2    NOF      Rb2SO4  UF4
BOH       CO              Cl2             H2O     NOF3     Rn      UF4+
BO2       CO+             Cl2O            H2O+    NO2      Rn+     UF4-
BO2-      COCl            Co              H2O2    NO2-     S       UF5
B_OH_2    COCl2           Co+             H2S     NO2Cl    S+      UF5+
BS        COFCl           Co-             H2SO4   NO2F     S-      UF5-
BS2       COF2            Cr              H2BOH   NO3      SCl     UF6
B2        COHCl           Cr+             HB_OH_2 NO3-     SCl2    UF6-
B2C       COHF            Cr-             H3BO3   NO3F     SCl2+   UO
B2Cl4     COS             CrN             H3B3O3  N2       SD      UO+
B2F4      CO2             CrO             H3B3O6  N2+      SF      UOF
B2H       CO2+            CrO2            H3F3    N2-      SF+     UOF2
B2H2      COOH            CrO3            H3O+    NCN      SF-     UOF3
B2H3      CP              CrO3-           H4F4    N2D2_cis SF2     UOF4
B2H3_db   CS              Cs              H5F5    N2F2     SF2+    UO2
B2H4      CS2             Cs+             H6F6    N2F4     SF2-    UO2+
B2H4_db   C2              Cs-             H7F7    N2H2     SF3     UO2-
B2H5      C2+             CsBO2           He      NH2NO2   SF3+    UO2F
B2H5_db   C2-             CsBr            He+     N2H4     SF3-    UO2F2
B2H6      C2Cl            CsCl            Hg      N2O      SF4     UO3
B2O       C2Cl2           CsF             Hg+     N2O+     SF4+    UO3-
B2O2      C2Cl3           CsH             HgBr2   N2O3     SF4-    V
B2O3      C2Cl4           CsI             I       N2O4     SF5     V+
B2_OH_4   C2Cl6           Csli            I+      N2O5     SF5+    V-
B2S       C2F             CsNO2           I-      N3       SF5-    VCl4
B2S2      C2FCl           CsNO3           IF5     N3H      SF6     VN
B2S3      C2FCl3          CsNa            IF7     Na       SF6-    VO
B3H7_C2v  C2F2            CsO             I2      Na+      SH      VO2
B3H7_Cs   C2F2Cl2         CsOH            In      Na-      SH-     V4O10
B3H9      C2F3            CsRb            In+     NaAlF4   SN      W
B3N3H6    C2F3Cl          Cs2             InBr    NaBO2    SO      W+
B3O3Cl3   C2F4            Cs2Br2          InBr2   NaBr     SO-     W-
B3O3FCl2  C2F6            Cs2CO3          InBr3   NaCN     SOF2    WCl6
B3O3F2Cl  C2H             Cs2Cl2          InCl    NaCl     SO2     WO
B3O3F3    C2HCl           Cs2F2           InCl2   NaF      SO2-    WOCl4
B4H4      C2HCl3          Cs2I2           InCl3   NaH      SO2Cl2  WO2
B4H10     C2HF            Cs2O            InF     NaI      SO2FCl  WO2Cl2
B4H12     C2HFCl2         Cs2O+           InF2    Nali     SO2F2   WO3
B5H9      C2HF2Cl         Cs2O2           InF3    NaNO2    SO3     WO3-
Ba        C2HF3           Cs2O2H2         InH     NaNO3    S2      Xe
Ba+       C2H2_vinylidene Cs2SO4          InI     NaO      S2-     Xe+
BaBr      C2H2Cl2         Cu              InI2    NaOH     S2Cl2   Zn
BaBr2     C2H2FCl         Cu+             InI3    NaOH+    S2F2    Zn+
BaCl      C2H2F2          Cu-             InO     Na2      S2O     Zr
BaCl+     CH2CO_ketene    CuCl            InOH    Na2Br2   S3      Zr+
BaCl2     O_CH_2O         CuF             In2Br2  Na2Cl2   S4      Zr-
BaF       HO_CO_2OH       CuF2            In2Br4  Na2F2    S5      ZrN
BaF+      C2H3_vinyl      CuO             In2Br6  Na2I2    S6      ZrO
BaF2      CH2Br-COOH      Cu2             In2Cl2  Na2O     S7      ZrO+
BaH       C2H3Cl          Cu3Cl3          In2Cl4  Na2O+    S8      ZrO2
BaI       CH2Cl-COOH      D               In2Cl6  Na2O2    Sc
BaI2      C2H3F           D+              In2F2   Na2O2H2  Sc+
BaO       CH3CN           D-              In2F4   Na2SO4   Sc-
BaO+      CH3CO_acetyl    DBr             In2F6   Na3Cl3   ScO
BaOH      C2H4            DCl             In2I2   Na3F3    ScO+
</pre></blockquote></html>"));
    end IdealGases;

    package Water "Medium models for water"
    extends Modelica.Icons.VariantsPackage;
    import Modelica.Media.Water.ConstantPropertyLiquidWater.simpleWaterConstants;

    package ConstantPropertyLiquidWater
      "Water: Simple liquid water medium (incompressible, constant data)"

      //   redeclare record extends FluidConstants
      //   end FluidConstants;

      constant Modelica.Media.Interfaces.Types.Basic.FluidConstants[1]
        simpleWaterConstants(
        each chemicalFormula="H2O",
        each structureFormula="H2O",
        each casRegistryNumber="7732-18-5",
        each iupacName="oxidane",
        each molarMass=0.018015268);

      extends Interfaces.PartialSimpleMedium(
        mediumName="SimpleLiquidWater",
        cp_const=4184,
        cv_const=4184,
        d_const=995.586,
        eta_const=1.e-3,
        lambda_const=0.598,
        a_const=1484,
        T_min=Cv.from_degC(-1),
        T_max=Cv.from_degC(130),
        T0=273.15,
        MM_const=0.018015268,
        fluidConstants=simpleWaterConstants);

      annotation (Documentation(info="<html>

</html>"));
    end ConstantPropertyLiquidWater;
    annotation (Documentation(info="<html>
<p>This package contains different medium models for water:</p>
<ul>
<li><strong>ConstantPropertyLiquidWater</strong><br>
    Simple liquid water medium (incompressible, constant data).</li>
<li><strong>IdealSteam</strong><br>
    Steam water medium as ideal gas from Media.IdealGases.SingleGases.H2O</li>
<li><strong>WaterIF97 derived models</strong><br>
    High precision water model according to the IAPWS/IF97 standard
    (liquid, steam, two phase region). Models with different independent
    variables are provided as well as models valid only
    for particular regions. The <strong>WaterIF97_ph</strong> model is valid
    in all regions and is the recommended one to use.</li>
</ul>
<h4>Overview of WaterIF97 derived water models</h4>
<p>
The WaterIF97 models calculate medium properties
for water in the <strong>liquid</strong>, <strong>gas</strong> and <strong>two phase</strong> regions
according to the IAPWS/IF97 standard, i.e., the accepted industrial standard
and best compromise between accuracy and computation time.
It has been part of the ThermoFluid Modelica library and been extended,
reorganized and documented to become part of the Modelica Standard library.</p>
<p>An important feature that distinguishes this implementation of the IF97 steam property standard
is that this implementation has been explicitly designed to work well in dynamic simulations. Computational
performance has been of high importance. This means that there often exist several ways to get the same result
from different functions if one of the functions is called often but can be optimized for that purpose.
</p>
<p>Three variable pairs can be the independent variables of the model:
</p>
<ol>
<li>Pressure <strong>p</strong> and specific enthalpy <strong>h</strong> are
    the most natural choice for general applications.
    This is the recommended choice for most general purpose
    applications, in particular for power plants.</li>
<li>Pressure <strong>p</strong> and temperature <strong>T</strong> are the most natural
    choice for applications where water is always in the same phase,
    both for liquid water and steam.</li>
<li>Density <strong>d</strong> and temperature <strong>T</strong> are explicit
    variables of the Helmholtz function in the near-critical
    region and can be the best choice for applications with
    super-critical or near-critical states.</li>
</ol>
<p>
The following quantities are always computed in Medium.BaseProperties:
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>Variable</strong></td>
      <td><strong>Unit</strong></td>
      <td><strong>Description</strong></td></tr>
  <tr><td>T</td>
      <td>K</td>
      <td>temperature</td></tr>
  <tr><td>u</td>
      <td>J/kg</td>
      <td>specific internal energy</td></tr>
  <tr><td>d</td>
      <td>kg/m^3</td>
      <td>density</td></tr>
  <tr><td>p</td>
      <td>Pa</td>
      <td>pressure</td></tr>
  <tr><td>h</td>
      <td>J/kg</td>
      <td>specific enthalpy</td></tr>
</table>
<p>
In some cases additional medium properties are needed.
A component that needs these optional properties has to call
one of the following functions:
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>Function call</strong></td>
      <td><strong>Unit</strong></td>
      <td><strong>Description</strong></td></tr>
  <tr><td>Medium.dynamicViscosity(medium.state)</td>
      <td>Pa.s</td>
      <td>dynamic viscosity</td></tr>
  <tr><td>Medium.thermalConductivity(medium.state)</td>
      <td>W/(m.K)</td>
      <td>thermal conductivity</td></tr>
  <tr><td>Medium.prandtlNumber(medium.state)</td>
      <td>1</td>
      <td>Prandtl number</td></tr>
  <tr><td>Medium.specificEntropy(medium.state)</td>
      <td>J/(kg.K)</td>
      <td>specific entropy</td></tr>
  <tr><td>Medium.heatCapacity_cp(medium.state)</td>
      <td>J/(kg.K)</td>
      <td>specific heat capacity at constant pressure</td></tr>
  <tr><td>Medium.heatCapacity_cv(medium.state)</td>
      <td>J/(kg.K)</td>
      <td>specific heat capacity at constant density</td></tr>
  <tr><td>Medium.isentropicExponent(medium.state)</td>
      <td>1</td>
      <td>isentropic exponent</td></tr>
  <tr><td>Medium.isentropicEnthalpy(pressure, medium.state)</td>
      <td>J/kg</td>
      <td>isentropic enthalpy</td></tr>
  <tr><td>Medium.velocityOfSound(medium.state)</td>
      <td>m/s</td>
      <td>velocity of sound</td></tr>
  <tr><td>Medium.isobaricExpansionCoefficient(medium.state)</td>
      <td>1/K</td>
      <td>isobaric expansion coefficient</td></tr>
  <tr><td>Medium.isothermalCompressibility(medium.state)</td>
      <td>1/Pa</td>
      <td>isothermal compressibility</td></tr>
  <tr><td>Medium.density_derp_h(medium.state)</td>
      <td>kg/(m3.Pa)</td>
      <td>derivative of density by pressure at constant enthalpy</td></tr>
  <tr><td>Medium.density_derh_p(medium.state)</td>
      <td>kg2/(m3.J)</td>
      <td>derivative of density by enthalpy at constant pressure</td></tr>
  <tr><td>Medium.density_derp_T(medium.state)</td>
      <td>kg/(m3.Pa)</td>
      <td>derivative of density by pressure at constant temperature</td></tr>
  <tr><td>Medium.density_derT_p(medium.state)</td>
      <td>kg/(m3.K)</td>
      <td>derivative of density by temperature at constant pressure</td></tr>
  <tr><td>Medium.density_derX(medium.state)</td>
      <td>kg/m3</td>
      <td>derivative of density by mass fraction</td></tr>
  <tr><td>Medium.molarMass(medium.state)</td>
      <td>kg/mol</td>
      <td>molar mass</td></tr>
</table>
<p>More details are given in
<a href=\"modelica://Modelica.Media.UsersGuide.MediumUsage.OptionalProperties\">
Modelica.Media.UsersGuide.MediumUsage.OptionalProperties</a>.

Many additional optional functions are defined to compute properties of
saturated media, either liquid (bubble point) or vapour (dew point).
The argument to such functions is a SaturationProperties record, which can be
set starting from either the saturation pressure or the saturation temperature.
With reference to a model defining a pressure p, a temperature T, and a
SaturationProperties record sat, the following functions are provided:
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><td><strong>Function call</strong></td>
      <td><strong>Unit</strong></td>
      <td><strong>Description</strong></td></tr>
  <tr><td>Medium.saturationPressure(T)</td>
      <td>Pa</td>
      <td>Saturation pressure at temperature T</td></tr>
  <tr><td>Medium.saturationTemperature(p)</td>
      <td>K</td>
      <td>Saturation temperature at pressure p</td></tr>
  <tr><td>Medium.saturationTemperature_derp(p)</td>
      <td>K/Pa</td>
      <td>Derivative of saturation temperature with respect to pressure</td></tr>
  <tr><td>Medium.bubbleEnthalpy(sat)</td>
      <td>J/kg</td>
      <td>Specific enthalpy at bubble point</td></tr>
  <tr><td>Medium.dewEnthalpy(sat)</td>
      <td>J/kg</td>
      <td>Specific enthalpy at dew point</td></tr>
  <tr><td>Medium.bubbleEntropy(sat)</td>
      <td>J/(kg.K)</td>
      <td>Specific entropy at bubble point</td></tr>
  <tr><td>Medium.dewEntropy(sat)</td>
      <td>J/(kg.K)</td>
      <td>Specific entropy at dew point</td></tr>
  <tr><td>Medium.bubbleDensity(sat)</td>
      <td>kg/m3</td>
      <td>Density at bubble point</td></tr>
  <tr><td>Medium.dewDensity(sat)</td>
      <td>kg/m3</td>
      <td>Density at dew point</td></tr>
  <tr><td>Medium.dBubbleDensity_dPressure(sat)</td>
      <td>kg/(m3.Pa)</td>
      <td>Derivative of density at bubble point with respect to pressure</td></tr>
  <tr><td>Medium.dDewDensity_dPressure(sat)</td>
      <td>kg/(m3.Pa)</td>
      <td>Derivative of density at dew point with respect to pressure</td></tr>
  <tr><td>Medium.dBubbleEnthalpy_dPressure(sat)</td>
      <td>J/(kg.Pa)</td>
      <td>Derivative of specific enthalpy at bubble point with respect to pressure</td></tr>
  <tr><td>Medium.dDewEnthalpy_dPressure(sat)</td>
      <td>J/(kg.Pa)</td>
      <td>Derivative of specific enthalpy at dew point with respect to pressure</td></tr>
  <tr><td>Medium.surfaceTension(sat)</td>
      <td>N/m</td>
      <td>Surface tension between liquid and vapour phase</td></tr>
</table>
<p>Details on usage and some examples are given in:
<a href=\"modelica://Modelica.Media.UsersGuide.MediumUsage.TwoPhase\">
Modelica.Media.UsersGuide.MediumUsage.TwoPhase</a>.
</p>
<p>Many further properties can be computed. Using the well-known Bridgman's Tables,
all first partial derivatives of the standard thermodynamic variables can be computed easily.
</p>
<p>
The documentation of the IAPWS/IF97 steam properties can be freely
distributed with computer implementations and are included here
(in directory Modelica/Resources/Documentation/Media/Water/IF97documentation):
</p>
<ul>
<li><a href=\"modelica://Modelica/Resources/Documentation/Media/Water/IF97documentation/IF97.pdf\">IF97.pdf</a> The standards document for the main part of the IF97.</li>
<li><a href=\"modelica://Modelica/Resources/Documentation/Media/Water/IF97documentation/Back3.pdf\">Back3.pdf</a> The backwards equations for region 3.</li>
<li><a href=\"modelica://Modelica/Resources/Documentation/Media/Water/IF97documentation/crits.pdf\">crits.pdf</a> The critical point data.</li>
<li><a href=\"modelica://Modelica/Resources/Documentation/Media/Water/IF97documentation/meltsub.pdf\">meltsub.pdf</a> The melting- and sublimation line formulation (not implemented)</li>
<li><a href=\"modelica://Modelica/Resources/Documentation/Media/Water/IF97documentation/surf.pdf\">surf.pdf</a> The surface tension standard definition</li>
<li><a href=\"modelica://Modelica/Resources/Documentation/Media/Water/IF97documentation/thcond.pdf\">thcond.pdf</a> The thermal conductivity standard definition</li>
<li><a href=\"modelica://Modelica/Resources/Documentation/Media/Water/IF97documentation/visc.pdf\">visc.pdf</a> The viscosity standard definition</li>
</ul>
</html>"));
    end Water;
  annotation (preferredView="info",Documentation(info="<html>
<p>
This library contains <a href=\"modelica://Modelica.Media.Interfaces\">interface</a>
definitions for media and the following <strong>property</strong> models for
single and multiple substance fluids with one and multiple phases:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Media.IdealGases\">Ideal gases:</a><br>
     1241 high precision gas models based on the
     NASA Glenn coefficients, plus ideal gas mixture models based
     on the same data.</li>
<li> <a href=\"modelica://Modelica.Media.Water\">Water models:</a><br>
     ConstantPropertyLiquidWater, WaterIF97 (high precision
     water model according to the IAPWS/IF97 standard)</li>
<li> <a href=\"modelica://Modelica.Media.Air\">Air models:</a><br>
     SimpleAir, DryAirNasa, ReferenceAir, MoistAir, ReferenceMoistAir.</li>
<li> <a href=\"modelica://Modelica.Media.Incompressible\">
     Incompressible media:</a><br>
     TableBased incompressible fluid models (properties are defined by tables rho(T),
     HeatCapacity_cp(T), etc.)</li>
<li> <a href=\"modelica://Modelica.Media.CompressibleLiquids\">
     Compressible liquids:</a><br>
     Simple liquid models with linear compressibility</li>
<li> <a href=\"modelica://Modelica.Media.R134a\">Refrigerant Tetrafluoroethane (R134a)</a>.</li>
</ul>
<p>
The following parts are useful, when newly starting with this library:</p>
<ul>
<li> <a href=\"modelica://Modelica.Media.UsersGuide\">Modelica.Media.UsersGuide</a>.</li>
<li> <a href=\"modelica://Modelica.Media.UsersGuide.MediumUsage\">Modelica.Media.UsersGuide.MediumUsage</a>
     describes how to use a medium model in a component model.</li>
<li> <a href=\"modelica://Modelica.Media.UsersGuide.MediumDefinition\">
     Modelica.Media.UsersGuide.MediumDefinition</a>
     describes how a new fluid medium model has to be implemented.</li>
<li> <a href=\"modelica://Modelica.Media.UsersGuide.ReleaseNotes\">Modelica.Media.UsersGuide.ReleaseNotes</a>
     summarizes the changes of the library releases.</li>
<li> <a href=\"modelica://Modelica.Media.Examples\">Modelica.Media.Examples</a>
     contains examples that demonstrate the usage of this library.</li>
</ul>
<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>February 01, 2017</em> by Thomas Beutlich:<br>
    Fixed data errors of the NASA Glenn coefficients in some ideal gases (CH2, CH3, CH3OOH, C2CL2, C2CL4, C2CL6, C2HCL, C2HCL3, CH2CO_ketene, O_CH_2O, HO_CO_2OH, CH2BrminusCOOH, C2H3CL, CH2CLminusCOOH, HO2, HO2minus, OD, ODminus), see <a href=\"https://github.com/modelica/ModelicaStandardLibrary/issues/1922\">#1922</a></li>
<li><em>May 16, 2013</em> by Stefan Wischhusen (XRG Simulation):<br>
    Added new media models Air.ReferenceMoistAir, Air.ReferenceAir, R134a.</li>
<li><em>May 25, 2011</em> by Francesco Casella:<br>Added min/max attributes to Water, TableBased, MixtureGasNasa, SimpleAir and MoistAir local types.</li>
<li><em>May 25, 2011</em> by Stefan Wischhusen:<br>Added individual settings for polynomial fittings of properties.</li>
</ul>
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

  package Thermal "Library of thermal system components to model heat transfer and simple thermo-fluid pipe flow"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;

    package HeatTransfer "Library of 1-dimensional heat transfer with lumped elements"
      extends Modelica.Icons.Package;

      package Components "Lumped thermal components"
        extends Modelica.Icons.Package;

        model HeatCapacitor "Lumped thermal element storing heat"
          parameter SI.HeatCapacity C
            "Heat capacity of element (= cp*m)";
          SI.Temperature T(start=293.15, displayUnit="degC")
            "Temperature of element";
          SI.TemperatureSlope der_T(start=0)
            "Time derivative of temperature (= der(T))";
          Interfaces.HeatPort_a port annotation (Placement(transformation(
                origin={0,-100},
                extent={{-10,-10},{10,10}},
                rotation=90)));
        equation
          T = port.T;
          der_T = der(T);
          C*der(T) = port.Q_flow;
          annotation (
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={
                Text(
                  extent={{-150,110},{150,70}},
                  textString="%name",
                  textColor={0,0,255}),
                Polygon(
                  points={{0,67},{-20,63},{-40,57},{-52,43},{-58,35},{-68,25},{-72,
                      13},{-76,-1},{-78,-15},{-76,-31},{-76,-43},{-76,-53},{-70,-65},
                      {-64,-73},{-48,-77},{-30,-83},{-18,-83},{-2,-85},{8,-89},{22,
                      -89},{32,-87},{42,-81},{54,-75},{56,-73},{66,-61},{68,-53},{
                      70,-51},{72,-35},{76,-21},{78,-13},{78,3},{74,15},{66,25},{54,
                      33},{44,41},{36,57},{26,65},{0,67}},
                  lineColor={160,160,164},
                  fillColor={192,192,192},
                  fillPattern=FillPattern.Solid),
                Polygon(
                  points={{-58,35},{-68,25},{-72,13},{-76,-1},{-78,-15},{-76,-31},{
                      -76,-43},{-76,-53},{-70,-65},{-64,-73},{-48,-77},{-30,-83},{-18,
                      -83},{-2,-85},{8,-89},{22,-89},{32,-87},{42,-81},{54,-75},{42,
                      -77},{40,-77},{30,-79},{20,-81},{18,-81},{10,-81},{2,-77},{-12,
                      -73},{-22,-73},{-30,-71},{-40,-65},{-50,-55},{-56,-43},{-58,-35},
                      {-58,-25},{-60,-13},{-60,-5},{-60,7},{-58,17},{-56,19},{-52,
                      27},{-48,35},{-44,45},{-40,57},{-58,35}},
                  fillColor={160,160,164},
                  fillPattern=FillPattern.Solid),
                Text(
                  extent={{-69,7},{71,-24}},
                  textString="%C")}),
            Documentation(info="<html>
<p>
This is a generic model for the heat capacity of a material.
No specific geometry is assumed beyond a total volume with
uniform temperature for the entire volume.
Furthermore, it is assumed that the heat capacity
is constant (independent of temperature).
</p>
<p>
The temperature T [Kelvin] of this component is a <strong>state</strong>.
A default of T = 25 degree Celsius (= Modelica.Units.Conversions.from_degC(25))
is used as start value for initialization.
This usually means that at start of integration the temperature of this
component is 25 degrees Celsius. You may, of course, define a different
temperature as start value for initialization. Alternatively, it is possible
to set parameter <strong>steadyStateStart</strong> to <strong>true</strong>. In this case
the additional equation '<strong>der</strong>(T) = 0' is used during
initialization, i.e., the temperature T is computed in such a way that
the component starts in <strong>steady state</strong>. This is useful in cases,
where one would like to start simulation in a suitable operating
point without being forced to integrate for a long time to arrive
at this point.
</p>
<p>
Note, that parameter <strong>steadyStateStart</strong> is not available in
the parameter menu of the simulation window, because its value
is utilized during translation to generate quite different
equations depending on its setting. Therefore, the value of this
parameter can only be changed before translating the model.
</p>
<p>
This component may be used for complicated geometries where
the heat capacity C is determined my measurements. If the component
consists mainly of one type of material, the <strong>mass m</strong> of the
component may be measured or calculated and multiplied with the
<strong>specific heat capacity cp</strong> of the component material to
compute C:
</p>
<blockquote><pre>
C = cp*m.
Typical values for cp at 20 degC in J/(kg.K):
   aluminium   896
   concrete    840
   copper      383
   iron        452
   silver      235
   steel       420 ... 500 (V2A)
   wood       2500
</pre></blockquote>
</html>"));
        end HeatCapacitor;

        model ThermalResistor
          "Lumped thermal element transporting heat without storing it"
          extends Interfaces.Element1D;
          parameter SI.ThermalResistance R
            "Constant thermal resistance of material";

        equation
          dT = R*Q_flow;
          annotation (
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={
                Rectangle(
                  extent={{-90,70},{90,-70}},
                  pattern=LinePattern.None,
                  fillColor={192,192,192},
                  fillPattern=FillPattern.Forward),
                Line(
                  points={{-90,70},{-90,-70}},
                  thickness=0.5),
                Line(
                  points={{90,70},{90,-70}},
                  thickness=0.5),
                Text(
                  extent={{-150,120},{150,78}},
                  textString="%name",
                  textColor={0,0,255}),
                Text(
                  extent={{-150,-80},{150,-110}},
                  textString="R=%R")}),
            Documentation(info="<html>
<p>
This is a model for transport of heat without storing it, same as the
<a href=\"modelica://Modelica.Thermal.HeatTransfer.Components.ThermalConductor\">ThermalConductor</a>
but using the thermal resistance instead of the thermal conductance as a parameter.
This is advantageous for series connections of ThermalResistors,
especially if it shall be allowed that a ThermalResistance is defined to be zero (i.e. no temperature difference).
</p>
</html>"));
        end ThermalResistor;

        model ThermalCollector "Collects m heat flows"
          parameter Integer m(min=1)=3 "Number of collected heat flows";
          Interfaces.HeatPort_a port_a[m]
            annotation (Placement(transformation(extent={{-10,110},{10,90}})));
          Interfaces.HeatPort_b port_b
            annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));

        equation
          port_b.Q_flow + sum(port_a.Q_flow) = 0;
          port_a.T = fill(port_b.T, m);
          annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={
                Text(
                  extent={{-150,-30},{150,-70}},
                  textString="%name",
                  textColor={0,0,255}),
                Text(
                  extent={{-150,80},{150,50}},
                  textString="m=%m"),
                Line(
                  points={{0,90},{0,40}},
                  color={181,0,0}),
                Rectangle(
                  extent={{-60,40},{60,30}},
                  lineColor={181,0,0},
                  fillColor={181,0,0},
                  fillPattern=FillPattern.Solid),
                Line(
                  points={{-60,30},{0,-30},{0,-90}},
                  color={181,0,0}),
                Line(
                  points={{0,-30},{-20,30}},
                  color={181,0,0}),
                Line(
                  points={{0,-30},{20,30}},
                  color={181,0,0}),
                Line(
                  points={{0,-30},{60,30}},
                  color={181,0,0})}),
            Documentation(info="<html>
<p>
This is a model to collect the heat flows from <em>m</em> heatports to one single heatport.
</p>
</html>"));
        end ThermalCollector;
        annotation (Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}), graphics={
          Rectangle(
            origin = {12,40},
            fillColor = {192,192,192},
            fillPattern = FillPattern.Backward,
            extent = {{-100,-100},{-70,18}}),
          Line(
            origin = {12,40},
            points = {{-44,16},{-44,-100}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{-4,16},{-4,-100}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{30,18},{30,-100}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{66,18},{66,-100}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{66,-100},{76,-80}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{66,-100},{56,-80}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{30,-100},{40,-80}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{30,-100},{20,-80}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{-4,-100},{6,-80}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{-4,-100},{-14,-80}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{-44,-100},{-34,-80}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{-44,-100},{-54,-80}},
            color = {0,127,255}),
          Line(
            origin = {12,40},
            points = {{-70,-60},{66,-60}},
            color = {191,0,0}),
          Line(
            origin = {12,40},
            points = {{46,-70},{66,-60}},
            color = {191,0,0}),
          Line(
            origin = {12,40},
            points = {{46,-50},{66,-60}},
            color = {191,0,0}),
          Line(
            origin = {12,40},
            points = {{46,-30},{66,-20}},
            color = {191,0,0}),
          Line(
            origin = {12,40},
            points = {{46,-10},{66,-20}},
            color = {191,0,0}),
          Line(
            origin = {12,40},
            points = {{-70,-20},{66,-20}},
            color = {191,0,0})}), Documentation(
              info="<html>

</html>"));
      end Components;

      package Sensors "Thermal sensors"
        extends Modelica.Icons.SensorsPackage;

        model TemperatureSensor "Absolute temperature sensor in Kelvin"

          Modelica.Blocks.Interfaces.RealOutput T(unit="K")
            "Absolute temperature as output signal"
            annotation (Placement(transformation(extent={{100,-10},{120,10}}), iconTransformation(extent={{100,-10},{120,10}})));
          Interfaces.HeatPort_a port annotation (Placement(transformation(extent={{
                    -110,-10},{-90,10}})));
        equation
          T = port.T;
          port.Q_flow = 0;
          annotation (
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={
                Ellipse(
                  extent={{-20,-98},{20,-60}},
                  lineThickness=0.5,
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-12,40},{12,-68}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid),
                Line(points={{12,0},{100,0}},color={0,0,127}),
                Line(points={{-90,0},{-12,0}}, color={191,0,0}),
                Polygon(
                  points={{-12,40},{-12,80},{-10,86},{-6,88},{0,90},{6,88},{10,86},
                      {12,80},{12,40},{-12,40}},
                  lineThickness=0.5),
                Line(
                  points={{-12,40},{-12,-64}},
                  thickness=0.5),
                Line(
                  points={{12,40},{12,-64}},
                  thickness=0.5),
                Line(points={{-40,-20},{-12,-20}}),
                Line(points={{-40,20},{-12,20}}),
                Line(points={{-40,60},{-12,60}}),
                Text(
                  extent={{-150,140},{150,100}},
                  textString="%name",
                  textColor={0,0,255}),
                Text(
                  extent={{20,60},{80,0}},
                  textColor={64,64,64},
                  textString="K")}),
            Documentation(info="<html>
<p>
This is an ideal absolute temperature sensor which returns
the temperature of the connected port in Kelvin as an output
signal.  The sensor itself has no thermal interaction with
whatever it is connected to.  Furthermore, no
thermocouple-like lags are associated with this
sensor model.
</p>
</html>"));
        end TemperatureSensor;

        model HeatFlowSensor "Heat flow rate sensor"
          extends Modelica.Icons.RoundSensor;
          Modelica.Blocks.Interfaces.RealOutput Q_flow(unit="W")
            "Heat flow from port_a to port_b as output signal" annotation (Placement(
                transformation(
                origin={0,-110},
                extent={{-10,-10},{10,10}},
                rotation=270), iconTransformation(
                extent={{-10,-10},{10,10}},
                rotation=270,
                origin={0,-110})));
          Interfaces.HeatPort_a port_a annotation (Placement(transformation(extent={{
                    -110,-10},{-90,10}})));
          Interfaces.HeatPort_b port_b annotation (Placement(transformation(extent={{
                    90,-10},{110,10}})));
        equation
          port_a.T = port_b.T;
          port_a.Q_flow + port_b.Q_flow = 0;
          Q_flow = port_a.Q_flow;
          annotation (
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={
                Line(points={{-70,0},{-90,0}}, color={191,0,0}),
                Line(points={{70,0},{90,0}}, color={191,0,0}),
                Line(points={{0,-70},{0,-100}},color={0,0,127}),
                Text(
                  extent={{-150,120},{150,80}},
                  textString="%name",
                  textColor={0,0,255}),
                Text(
                  extent={{-30,-10},{30,-70}},
                  textColor={64,64,64},
                  textString="W")}),
            Documentation(info="<html>
<p>
This model is capable of monitoring the heat flow rate flowing through
this component. The sensed value of heat flow rate is the amount that
passes through this sensor while keeping the temperature drop across the
sensor zero.  This is an ideal model so it does not absorb any energy
and it has no direct effect on the thermal response of a system it is included in.
The output signal is positive, if the heat flows from port_a to port_b.
</p>
</html>"));
        end HeatFlowSensor;
        annotation (Documentation(info="<html>

</html>"));
      end Sensors;

      package Sources "Thermal sources"
        extends Modelica.Icons.SourcesPackage;

        model FixedTemperature "Fixed temperature boundary condition in Kelvin"

          parameter SI.Temperature T "Fixed temperature at port";
          Interfaces.HeatPort_b port annotation (Placement(transformation(extent={{90,
                    -10},{110,10}})));
        equation
          port.T = T;
          annotation (
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={
                Text(
                  extent={{-150,150},{150,110}},
                  textString="%name",
                  textColor={0,0,255}),
                Text(
                  extent={{-150,-110},{150,-140}},
                  textString="T=%T"),
                Rectangle(
                  extent={{-100,100},{100,-100}},
                  pattern=LinePattern.None,
                  fillColor={159,159,223},
                  fillPattern=FillPattern.Backward),
                Text(
                  extent={{0,0},{-100,-100}},
                  textString="K"),
                Line(
                  points={{-52,0},{56,0}},
                  color={191,0,0},
                  thickness=0.5),
                Polygon(
                  points={{50,-20},{50,20},{90,0},{50,-20}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid)}),
            Documentation(info="<html>
<p>
This model defines a fixed temperature T at its port in Kelvin,
i.e., it defines a fixed temperature as a boundary condition.
</p>
</html>"));
        end FixedTemperature;

        model PrescribedTemperature
          "Variable temperature boundary condition in Kelvin"

          Interfaces.HeatPort_b port annotation (Placement(transformation(extent={{90,
                    -10},{110,10}})));
          Modelica.Blocks.Interfaces.RealInput T(unit="K") annotation (Placement(transformation(
                  extent={{-140,-20},{-100,20}})));
        equation
          port.T = T;
          annotation (
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={
                Rectangle(
                  extent={{-100,100},{100,-100}},
                  pattern=LinePattern.None,
                  fillColor={159,159,223},
                  fillPattern=FillPattern.Backward),
                Line(
                  points={{-102,0},{64,0}},
                  color={191,0,0},
                  thickness=0.5),
                Text(
                  extent={{0,0},{-100,-100}},
                  textString="K"),
                Text(
                  extent={{-150,150},{150,110}},
                  textString="%name",
                  textColor={0,0,255}),
                Polygon(
                  points={{50,-20},{50,20},{90,0},{50,-20}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid)}),
            Documentation(info="<html>
<p>
This model represents a variable temperature boundary condition.
The temperature in [K] is given as input signal <strong>T</strong>
to the model. The effect is that an instance of this model acts as
an infinite reservoir able to absorb or generate as much energy
as required to keep the temperature at the specified value.
</p>
</html>"));
        end PrescribedTemperature;

        model PrescribedHeatFlow "Prescribed heat flow boundary condition"
          parameter SI.Temperature T_ref=293.15
            "Reference temperature";
          parameter SI.LinearTemperatureCoefficient alpha=0
            "Temperature coefficient of heat flow rate";
          Modelica.Blocks.Interfaces.RealInput Q_flow(unit="W")
                annotation (Placement(transformation(
                origin={-100,0},
                extent={{20,-20},{-20,20}},
                rotation=180)));
          Interfaces.HeatPort_b port annotation (Placement(transformation(extent={{90,
                    -10},{110,10}})));
        equation
          port.Q_flow = -Q_flow*(1 + alpha*(port.T - T_ref));
          annotation (
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={
                Line(
                  points={{-60,-20},{40,-20}},
                  color={191,0,0},
                  thickness=0.5),
                Line(
                  points={{-60,20},{40,20}},
                  color={191,0,0},
                  thickness=0.5),
                Line(
                  points={{-80,0},{-60,-20}},
                  color={191,0,0},
                  thickness=0.5),
                Line(
                  points={{-80,0},{-60,20}},
                  color={191,0,0},
                  thickness=0.5),
                Polygon(
                  points={{40,0},{40,40},{70,20},{40,0}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid),
                Polygon(
                  points={{40,-40},{40,0},{70,-20},{40,-40}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{70,40},{90,-40}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid),
                Text(
                  extent={{-150,100},{150,60}},
                  textString="%name",
                  textColor={0,0,255})}),
            Documentation(info="<html>
<p>
This model allows a specified amount of heat flow rate to be \"injected\"
into a thermal system at a given port.  The amount of heat
is given by the input signal Q_flow into the model. The heat flows into the
component to which the component PrescribedHeatFlow is connected,
if the input signal is positive.
</p>
<p>
If parameter alpha is &lt;&gt; 0, the heat flow is multiplied by (1 + alpha*(port.T - T_ref))
in order to simulate temperature dependent losses (which are given with respect to reference temperature T_ref).
</p>
</html>"));
        end PrescribedHeatFlow;
      end Sources;

      package Interfaces "Connectors and partial models"
        extends Modelica.Icons.InterfacesPackage;

        partial connector HeatPort "Thermal port for 1-dim. heat transfer"
          SI.Temperature T "Port temperature";
          flow SI.HeatFlowRate Q_flow
            "Heat flow rate (positive if flowing from outside into the component)";
          annotation (Documentation(info="<html>

</html>"));
        end HeatPort;

        connector HeatPort_a
          "Thermal port for 1-dim. heat transfer (filled rectangular icon)"

          extends HeatPort;

          annotation(defaultComponentName = "port_a",
            Documentation(info="<html>
<p>This connector is used for 1-dimensional heat flow between components.
The variables in the connector are:</p>
<blockquote><pre>
T       Temperature in [Kelvin].
Q_flow  Heat flow rate in [Watt].
</pre></blockquote>
<p>According to the Modelica sign convention, a <strong>positive</strong> heat flow
rate <strong>Q_flow</strong> is considered to flow <strong>into</strong> a component. This
convention has to be used whenever this connector is used in a model
class.</p>
<p>Note, that the two connector classes <strong>HeatPort_a</strong> and
<strong>HeatPort_b</strong> are identical with the only exception of the different
<strong>icon layout</strong>.</p></html>"),     Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid)}),
            Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}), graphics={Rectangle(
                  extent={{-50,50},{50,-50}},
                  lineColor={191,0,0},
                  fillColor={191,0,0},
                  fillPattern=FillPattern.Solid), Text(
                  extent={{-120,120},{100,60}},
                  textColor={191,0,0},
                  textString="%name")}));
        end HeatPort_a;

        connector HeatPort_b
          "Thermal port for 1-dim. heat transfer (unfilled rectangular icon)"

          extends HeatPort;

          annotation(defaultComponentName = "port_b",
            Documentation(info="<html>
<p>This connector is used for 1-dimensional heat flow between components.
The variables in the connector are:</p>
<blockquote><pre>
T       Temperature in [Kelvin].
Q_flow  Heat flow rate in [Watt].
</pre></blockquote>
<p>According to the Modelica sign convention, a <strong>positive</strong> heat flow
rate <strong>Q_flow</strong> is considered to flow <strong>into</strong> a component. This
convention has to be used whenever this connector is used in a model
class.</p>
<p>Note, that the two connector classes <strong>HeatPort_a</strong> and
<strong>HeatPort_b</strong> are identical with the only exception of the different
<strong>icon layout</strong>.</p></html>"),     Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
                    {100,100}}), graphics={Rectangle(
                  extent={{-50,50},{50,-50}},
                  lineColor={191,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid), Text(
                  extent={{-100,120},{120,60}},
                  textColor={191,0,0},
                  textString="%name")}),
            Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
                    100,100}}), graphics={Rectangle(
                  extent={{-100,100},{100,-100}},
                  lineColor={191,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid)}));
        end HeatPort_b;

        partial model Element1D
          "Partial heat transfer element with two HeatPort connectors that does not store energy"

          SI.HeatFlowRate Q_flow
            "Heat flow rate from port_a -> port_b";
          SI.TemperatureDifference dT "port_a.T - port_b.T";
      public
          HeatPort_a port_a annotation (Placement(transformation(extent={{-110,-10},
                    {-90,10}})));
          HeatPort_b port_b annotation (Placement(transformation(extent={{90,-10},{
                    110,10}})));
        equation
          dT = port_a.T - port_b.T;
          port_a.Q_flow = Q_flow;
          port_b.Q_flow = -Q_flow;
          annotation (Documentation(info="<html>
<p>
This partial model contains the basic connectors and variables to
allow heat transfer models to be created that <strong>do not store energy</strong>,
This model defines and includes equations for the temperature
drop across the element, <strong>dT</strong>, and the heat flow rate
through the element from port_a to port_b, <strong>Q_flow</strong>.
</p>
<p>
By extending this model, it is possible to write simple
constitutive equations for many types of heat transfer components.
</p>
</html>"));
        end Element1D;
        annotation (Documentation(info="<html>

</html>"));
      end Interfaces;
      annotation (
         Icon(coordinateSystem(preserveAspectRatio = true, extent = {{-100,-100},{100,100}}), graphics={
          Polygon(
            origin = {13.758,27.517},
            lineColor = {128,128,128},
            fillColor = {192,192,192},
            fillPattern = FillPattern.Solid,
            points = {{-54,-6},{-61,-7},{-75,-15},{-79,-24},{-80,-34},{-78,-42},{-73,-49},{-64,-51},{-57,-51},{-47,-50},{-41,-43},{-38,-35},{-40,-27},{-40,-20},{-42,-13},{-47,-7},{-54,-5},{-54,-6}}),
        Polygon(
            origin = {13.758,27.517},
            fillColor = {160,160,164},
            fillPattern = FillPattern.Solid,
            points = {{-75,-15},{-79,-25},{-80,-34},{-78,-42},{-72,-49},{-64,-51},{-57,-51},{-47,-50},{-57,-47},{-65,-45},{-71,-40},{-74,-33},{-76,-23},{-75,-15},{-75,-15}}),
          Polygon(
            origin = {13.758,27.517},
            lineColor = {160,160,164},
            fillColor = {192,192,192},
            fillPattern = FillPattern.Solid,
            points = {{39,-6},{32,-7},{18,-15},{14,-24},{13,-34},{15,-42},{20,-49},{29,-51},{36,-51},{46,-50},{52,-43},{55,-35},{53,-27},{53,-20},{51,-13},{46,-7},{39,-5},{39,-6}}),
          Polygon(
            origin = {13.758,27.517},
            fillColor = {160,160,164},
            fillPattern = FillPattern.Solid,
            points = {{18,-15},{14,-25},{13,-34},{15,-42},{21,-49},{29,-51},{36,-51},{46,-50},{36,-47},{28,-45},{22,-40},{19,-33},{17,-23},{18,-15},{18,-15}}),
          Polygon(
            origin = {13.758,27.517},
            lineColor = {191,0,0},
            fillColor = {191,0,0},
            fillPattern = FillPattern.Solid,
            points = {{-9,-23},{-9,-10},{18,-17},{-9,-23}}),
          Line(
            origin = {13.758,27.517},
            points = {{-41,-17},{-9,-17}},
            color = {191,0,0},
            thickness = 0.5),
          Line(
            origin = {13.758,27.517},
            points = {{-17,-40},{15,-40}},
            color = {191,0,0},
            thickness = 0.5),
          Polygon(
            origin = {13.758,27.517},
            lineColor = {191,0,0},
            fillColor = {191,0,0},
            fillPattern = FillPattern.Solid,
            points = {{-17,-46},{-17,-34},{-40,-40},{-17,-46}})}),
                                Documentation(info="<html>
<p>
This package contains components to model <strong>1-dimensional heat transfer</strong>
with lumped elements.</p>
</html>",     revisions="<html>
<ul>
<li><em>July 15, 2002</em>
       by Michael Tiller, <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Nikolaus Sch&uuml;rmann:<br>
       Implemented.
</li>
<li><em>June 13, 2005</em>
       by <a href=\"https://www.haumer.at/\">Anton Haumer</a><br>
       Refined placing of connectors (cosmetic).<br>
       Refined all Examples; removed Examples.FrequencyInverter, introducing Examples.Motor<br>
       Introduced temperature dependent correction (1 + alpha*(T - T_ref)) in Fixed/PrescribedHeatFlow<br>
</li>
  <li> v1.1.1 2007/11/13 Anton Haumer<br>
       components moved to sub-packages</li>
  <li> v1.2.0 2009/08/26 Anton Haumer<br>
       added component ThermalCollector</li>

</ul>
</html>"));
    end HeatTransfer;
    annotation (
     Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
      Line(
      origin={-47.5,11.6667},
      points={{-2.5,-91.6667},{17.5,-71.6667},{-22.5,-51.6667},{17.5,-31.6667},{-22.5,-11.667},{17.5,8.3333},{-2.5,28.3333},{-2.5,48.3333}},
        smooth=Smooth.Bezier),
      Polygon(
      origin={-50.0,68.333},
      pattern=LinePattern.None,
      fillPattern=FillPattern.Solid,
        points={{0.0,21.667},{-10.0,-8.333},{10.0,-8.333}}),
      Line(
      origin={2.5,11.6667},
      points={{-2.5,-91.6667},{17.5,-71.6667},{-22.5,-51.6667},{17.5,-31.6667},{-22.5,-11.667},{17.5,8.3333},{-2.5,28.3333},{-2.5,48.3333}},
        smooth=Smooth.Bezier),
      Polygon(
      origin={0.0,68.333},
      pattern=LinePattern.None,
      fillPattern=FillPattern.Solid,
        points={{0.0,21.667},{-10.0,-8.333},{10.0,-8.333}}),
      Line(
      origin={52.5,11.6667},
      points={{-2.5,-91.6667},{17.5,-71.6667},{-22.5,-51.6667},{17.5,-31.6667},{-22.5,-11.667},{17.5,8.3333},{-2.5,28.3333},{-2.5,48.3333}},
        smooth=Smooth.Bezier),
      Polygon(
      origin={50.0,68.333},
      pattern=LinePattern.None,
      fillPattern=FillPattern.Solid,
        points={{0.0,21.667},{-10.0,-8.333},{10.0,-8.333}})}),
      Documentation(info="<html>
<p>
This package contains libraries to model heat transfer
and fluid heat flow.
</p>
</html>"));
  end Thermal;

  package Math "Library of mathematical functions (e.g., sin, cos) and of functions operating on vectors and matrices"
    extends Modelica.Icons.Package;

    package Nonlinear "Library of functions operating on nonlinear equations"
      extends Modelica.Icons.Package;

      package Interfaces "Interfaces for functions"
        extends Modelica.Icons.InterfacesPackage;

      encapsulated partial function partialScalarFunction
          "Interface for a function with one input and one output Real signal"
        import Modelica;
        extends Modelica.Icons.Function;
        input Real u "Independent variable";
        output Real y "Dependent variable y=f(u)";
          annotation (Documentation(info="<html>
<p>
This partial function defines the interface of a function with
one input and one output Real signal. The scalar functions
of <a href=\"modelica://Modelica.Math.Nonlinear\">Modelica.Math.Nonlinear</a>
are derived from this base type by inheritance.
This allows to use these functions directly as function arguments
to a function, see, .e.g.,
<a href=\"modelica://Modelica.Math.Nonlinear.Examples\">Math.Nonlinear.Examples</a>.
</p>

</html>"));
      end partialScalarFunction;
        annotation (Documentation(info="<html>
<p>
Interface definitions of functions. The main purpose is to use functions
derived from these interface definitions as function arguments
to a function, see, .e.g.,
<a href=\"modelica://Modelica.Math.Nonlinear.Examples\">Math.Nonlinear.Examples</a>.
</p>
</html>"));
      end Interfaces;

      function solveOneNonlinearEquation
        "Solve f(u) = 0 in a very reliable and efficient way (f(u_min) and f(u_max) must have different signs)"
        extends Modelica.Icons.Function;
        import Modelica.Utilities.Streams.error;

        input Modelica.Math.Nonlinear.Interfaces.partialScalarFunction f
          "Function y = f(u); u is computed so that y=0";
        input Real u_min "Lower bound of search interval";
        input Real u_max "Upper bound of search interval";
        input Real tolerance=100*Modelica.Constants.eps
          "Relative tolerance of solution u";
        output Real u "Value of independent variable u so that f(u) = 0";

    protected
        constant Real eps=Modelica.Constants.eps "Machine epsilon";
        Real a=u_min "Current best minimum interval value";
        Real b=u_max "Current best maximum interval value";
        Real c "Intermediate point a <= c <= b";
        Real d;
        Real e "b - a";
        Real m;
        Real s;
        Real p;
        Real q;
        Real r;
        Real tol;
        Real fa "= f(a)";
        Real fb "= f(b)";
        Real fc;
        Boolean found=false;
      algorithm
        // Check that f(u_min) and f(u_max) have different sign
        fa := f(u_min);
        fb := f(u_max);
        fc := fb;
        if fa > 0.0 and fb > 0.0 or fa < 0.0 and fb < 0.0 then
          error(
            "The arguments u_min and u_max provided in the function call\n"+
            "    solveOneNonlinearEquation(f,u_min,u_max)\n" +
            "do not bracket the root of the single non-linear equation 0=f(u):\n" +
            "  u_min  = " + String(u_min) + "\n" + "  u_max  = " + String(u_max)
             + "\n" + "  fa = f(u_min) = " + String(fa) + "\n" +
            "  fb = f(u_max) = " + String(fb) + "\n" +
            "fa and fb must have opposite sign which is not the case");
        end if;

        // Initialize variables
        c := a;
        fc := fa;
        e := b - a;
        d := e;

        // Search loop
        while not found loop
          if abs(fc) < abs(fb) then
            a := b;
            b := c;
            c := a;
            fa := fb;
            fb := fc;
            fc := fa;
          end if;

          tol := 2*eps*abs(b) + tolerance;
          m := (c - b)/2;

          if abs(m) <= tol or fb == 0.0 then
            // root found (interval is small enough)
            found := true;
            u := b;
          else
            // Determine if a bisection is needed
            if abs(e) < tol or abs(fa) <= abs(fb) then
              e := m;
              d := e;
            else
              s := fb/fa;
              if a == c then
                // linear interpolation
                p := 2*m*s;
                q := 1 - s;
              else
                // inverse quadratic interpolation
                q := fa/fc;
                r := fb/fc;
                p := s*(2*m*q*(q - r) - (b - a)*(r - 1));
                q := (q - 1)*(r - 1)*(s - 1);
              end if;

              if p > 0 then
                q := -q;
              else
                p := -p;
              end if;

              s := e;
              e := d;
              if 2*p < 3*m*q - abs(tol*q) and p < abs(0.5*s*q) then
                // interpolation successful
                d := p/q;
              else
                // use bi-section
                e := m;
                d := e;
              end if;
            end if;

            // Best guess value is defined as "a"
            a := b;
            fa := fb;
            b := b + (if abs(d) > tol then d else if m > 0 then tol else -tol);
            fb := f(b);

            if fb > 0 and fc > 0 or fb < 0 and fc < 0 then
              // initialize variables
              c := a;
              fc := fa;
              e := b - a;
              d := e;
            end if;
          end if;
        end while;

        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
<strong>solveOneNonlinearEquation</strong>(function f(), u_min, u_max);
<strong>solveOneNonlinearEquation</strong>(function f(), u_min, u_max, tolerance=100*Modelica.Constants.eps);
</pre></blockquote>

<h4>Description</h4>

<p>
This function determines the solution of <strong>one non-linear algebraic equation</strong> \"y=f(u)\"
in <strong>one unknown</strong> \"u\" in a reliable way. It is one of the best numerical
algorithms for this purpose. As input, the nonlinear function f(u)
has to be given, as well as an interval u_min, u_max that
contains the solution, i.e., \"f(u_min)\" and \"f(u_max)\" must
have a different sign. The function computes a smaller interval
in which a sign change is present using the relative tolerance
\"tolerance\" that can be given as 4th input argument.
</p>

<p>
The interval reduction is performed using
inverse quadratic interpolation (interpolating with a quadratic polynomial
through the last 3 points and computing the zero). If this fails,
bisection is used, which always reduces the interval by a factor of 2.
The inverse quadratic interpolation method has superlinear convergence.
This is roughly the same convergence rate as a globally convergent Newton
method, but without the need to compute derivatives of the non-linear
function. The solver function is a direct mapping of the Algol 60 procedure
\"zero\" to Modelica, from:
</p>

<blockquote>
<dl>
<dt> Brent R.P.:</dt>
<dd> <strong>Algorithms for Minimization without derivatives</strong>.
     Prentice Hall, 1973, pp. 58-59.<br>
     Download: <a href=\"https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf\">https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf</a><br>
     Errata and new print: <a href=\"https://maths-people.anu.edu.au/~brent/pub/pub011.html\">https://maths-people.anu.edu.au/~brent/pub/pub011.html</a>
</dd>
</dl>
</blockquote>

<h4>Example</h4>

<p>
See the examples in <a href=\"modelica://Modelica.Math.Nonlinear.Examples\">Modelica.Math.Nonlinear.Examples</a>.
</p>
</html>"));
      end solveOneNonlinearEquation;
      annotation (Documentation(info="<html>
<p>
This package contains functions to perform tasks such as numerically integrating
a function, or solving a nonlinear algebraic equation system.
The common feature of the functions in this package is
that the nonlinear characteristics are passed as user definable
functions.
</p>

<p>
For details about how to define and to use functions as input arguments
to functions, see
<a href=\"modelica://ModelicaReference.Classes.'function'\">ModelicaReference.Classes.'function'</a>
or <a href=\"https://specification.modelica.org/v3.4/Ch12.html#functional-input-arguments-to-functions\">Section 12.4.2
(Functional Input Arguments to Functions) of the Modelica 3.4 specification</a>.
</p>

</html>",     revisions="<html>
<ul>
<li><em>July 2010 </em> by Martin Otter (DLR-RM):<br>
    Included in MSL 3.2, adapted, and documentation improved</li>

<li><em>March 2010 </em> by Andreas Pfeiffer (DLR-RM):<br>
    Adapted the quadrature function from Gerhard Schillhuber and
    the solution of one non-linear equation in one unknown from
    Modelica.Media.Common.OneNonLinearEquation so that
    function objects are used.</li>

<li><em>June 2002 </em> by Gerhard Schillhuber (master thesis at DLR-RM):<br>
       Adaptive quadrature to compute the curve length of a Spline.</li>
</ul>
</html>"),     Icon(graphics={Polygon(points={{-44,-52},{-44,-26},{-17.1,
                  44.4},{-11.4,52.6},{-5.8,57.1},{-0.2,57.8},{5.4,54.6},{11.1,47.7},
                  {16.7,37.4},{23.1,22.1},{31.17,-0.8},{48,-52},{-44,-52}},
              lineColor={135,135,135},
              fillColor={215,215,215},
              fillPattern=FillPattern.Solid)}));
    end Nonlinear;

        package Polynomials "Library of functions operating on polynomials (including polynomial fitting)"
          extends Modelica.Icons.FunctionsPackage;

          function evaluate "Evaluate polynomial at a given abscissa value"
            extends Modelica.Icons.Function;
            input Real p[:]
              "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real u "Abscissa value";
            output Real y "Value of polynomial at u";
          algorithm
            y := p[1];
            for j in 2:size(p, 1) loop
              y := p[j] + u*y;
            end for;
            annotation(derivative(zeroDerivative=p)=evaluate_der);
          end evaluate;

          function evaluateWithRange
            "Evaluate polynomial at a given abscissa value with linear extrapolation outside of the defined range"
            extends Modelica.Icons.Function;
            input Real p[:]
              "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real uMin "Polynomial valid in the range uMin .. uMax";
            input Real uMax "Polynomial valid in the range uMin .. uMax";
            input Real u "Abscissa value";
            output Real y
              "Value of polynomial at u. Outside of uMin,uMax, linear extrapolation is used";
          algorithm
            if u < uMin then
              y := evaluate(p, uMin) - evaluate_der(
                      p,
                      uMin,
                      uMin - u);
            elseif u > uMax then
              y := evaluate(p, uMax) + evaluate_der(
                      p,
                      uMax,
                      u - uMax);
            else
              y := evaluate(p, u);
            end if;
            annotation (derivative(
                zeroDerivative=p,
                zeroDerivative=uMin,
                zeroDerivative=uMax) = evaluateWithRange_der);
          end evaluateWithRange;

          function evaluate_der
            "Evaluate derivative of polynomial at a given abscissa value"
            extends Modelica.Icons.Function;
            input Real p[:]
              "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real u "Abscissa value";
            input Real du "Delta of abscissa value";
            output Real dy "Value of derivative of polynomial at u";
        protected
            Integer n=size(p, 1);
          algorithm
            dy := p[1]*(n - 1);
            for j in 2:size(p, 1)-1 loop
              dy := p[j]*(n - j) + u*dy;
            end for;
            dy := dy*du;
          end evaluate_der;

          function evaluateWithRange_der
            "Evaluate derivative of polynomial at a given abscissa value with extrapolation outside of the defined range"
            extends Modelica.Icons.Function;
            input Real p[:]
              "Polynomial coefficients (p[1] is coefficient of highest power)";
            input Real uMin "Polynomial valid in the range uMin .. uMax";
            input Real uMax "Polynomial valid in the range uMin .. uMax";
            input Real u "Abscissa value";
            input Real du "Delta of abscissa value";
            output Real dy "Value of derivative of polynomial at u";
          algorithm
            if u < uMin then
              dy := evaluate_der(
                      p,
                      uMin,
                      du);
            elseif u > uMax then
              dy := evaluate_der(
                      p,
                      uMax,
                      du);
            else
              dy := evaluate_der(
                      p,
                      u,
                      du);
            end if;
          end evaluateWithRange_der;
          annotation (Documentation(info="<html>
<p>
This package contains functions to operate on polynomials,
in particular to determine the derivative and the integral
of a polynomial and to use a polynomial to fit a given set
of data points.
</p>

<p>
Copyright &copy; 2004-2020, Modelica Association and contributors
</p>
</html>",         revisions="<html>
<ul>
<li><em>Oct. 22, 2004</em> by Martin Otter (DLR):<br>
       Renamed functions to not have abbreviations.<br>
       Based fitting on LAPACK<br>
       New function to return the polynomial of an indefinite integral</li>
<li><em>Sept. 3, 2004</em> by Jonas Eborn (Scynamics):<br>
       polyderval, polyintval added</li>
<li><em>March 1, 2004</em> by Martin Otter (DLR):<br>
       first version implemented</li>
</ul>
</html>"));
        end Polynomials;

  package Icons "Icons for Math"
    extends Modelica.Icons.IconsPackage;

    partial function AxisLeft
      "Basic icon for mathematical function with y-axis on left side"

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                100}}), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{-80,-80},{-80,68}}, color={192,192,192}),
            Polygon(
              points={{-80,90},{-88,68},{-72,68},{-80,90}},
              lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              textColor={0,0,255})}),
        Documentation(info="<html>
<p>
Icon for a mathematical function, consisting of an y-axis on the left side.
It is expected, that an x-axis is added and a plot of the function.
</p>
</html>"));
    end AxisLeft;

    partial function AxisCenter
      "Basic icon for mathematical function with y-axis in the center"

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
                100}}), graphics={
            Rectangle(
              extent={{-100,100},{100,-100}},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Line(points={{0,-80},{0,68}}, color={192,192,192}),
            Polygon(
              points={{0,90},{-8,68},{8,68},{0,90}},
              lineColor={192,192,192},
              fillColor={192,192,192},
              fillPattern=FillPattern.Solid),
            Text(
              extent={{-150,150},{150,110}},
              textString="%name",
              textColor={0,0,255})}),
        Documentation(info="<html>
<p>
Icon for a mathematical function, consisting of an y-axis in the middle.
It is expected, that an x-axis is added and a plot of the function.
</p>
</html>"));
    end AxisCenter;
  end Icons;

  function isEqual "Determine if two Real scalars are numerically identical"
    extends Modelica.Icons.Function;
    input Real s1 "First scalar";
    input Real s2 "Second scalar";
    input Real eps(min=0) = 0
      "The two scalars are identical if abs(s1-s2) <= eps";
    output Boolean result "= true, if scalars are identical";
  algorithm
    result := abs(s1 - s2) <= eps;
    annotation (Inline=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Math.<strong>isEqual</strong>(s1, s2);
Math.<strong>isEqual</strong>(s1, s2, eps=0);
</pre></blockquote>
<h4>Description</h4>
<p>
The function call \"<code>Math.isEqual(s1, s2)</code>\" returns <strong>true</strong>,
if the two Real scalars s1 and s2 are identical. Otherwise the function
returns <strong>false</strong>. The equality check is performed by
\"abs(s1-s2) &le; eps\", where \"eps\"
can be provided as third argument of the function. Default is \"eps = 0\".
</p>
<h4>Example</h4>
<blockquote><pre>
  Real s1 = 2.0;
  Real s2 = 2.0;
  Real s3 = 2.000001;
  Boolean result;
<strong>algorithm</strong>
  result := Math.isEqual(s1,s2);     // = <strong>true</strong>
  result := Math.isEqual(s1,s3);     // = <strong>false</strong>
  result := Math.isEqual(s1,s3,0.1); // = <strong>true</strong>
</pre></blockquote>
<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Math.Vectors.isEqual\">Vectors.isEqual</a>,
<a href=\"modelica://Modelica.Math.Matrices.isEqual\">Matrices.isEqual</a>,
<a href=\"modelica://Modelica.Utilities.Strings.isEqual\">Strings.isEqual</a>
</p>
</html>"));
  end isEqual;

  function cos "Cosine"
    extends Modelica.Math.Icons.AxisLeft;
    input Modelica.Units.SI.Angle u "Independent variable";
    output Real y "Dependent variable y=cos(u)";

  external "builtin" y = cos(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,0},{68,0}}, color={192,192,192}),
          Polygon(
            points={{90,0},{68,8},{68,-8},{90,0}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,80},{-74.4,78.1},{-68.7,72.3},{-63.1,63},{-56.7,48.7},
                {-48.6,26.6},{-29.3,-32.5},{-22.1,-51.7},{-15.7,-65.3},{-10.1,-73.8},
                {-4.42,-78.8},{1.21,-79.9},{6.83,-77.1},{12.5,-70.6},{18.1,-60.6},
                {24.5,-45.7},{32.6,-23},{50.3,31.3},{57.5,50.7},{63.9,64.6},{69.5,
                73.4},{75.2,78.6},{80,80}}),
          Text(
            extent={{-36,82},{36,34}},
            textColor={192,192,192},
            textString="cos")}),
      Documentation(info="<html>
<p>
This function returns y = cos(u), with -&infin; &lt; u &lt; &infin;:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/cos.png\">
</p>
</html>"));
  end cos;

  function tan "Tangent (u shall not be -pi/2, pi/2, 3*pi/2, ...)"
    extends Modelica.Math.Icons.AxisCenter;
    input Modelica.Units.SI.Angle u "Independent variable";
    output Real y "Dependent variable y=tan(u)";

  external "builtin" y = tan(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,0},{68,0}}, color={192,192,192}),
          Polygon(
            points={{90,0},{68,8},{68,-8},{90,0}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-80},{-78.4,-68.4},{-76.8,-59.7},{-74.4,-50},{-71.2,-40.9},
                {-67.1,-33},{-60.7,-24.8},{-51.1,-17.2},{-35.8,-9.98},{-4.42,-1.07},
                {33.4,9.12},{49.4,16.2},{59.1,23.2},{65.5,30.6},{70.4,39.1},{73.6,
                47.4},{76,56.1},{77.6,63.8},{80,80}}),
          Text(
            extent={{-90,72},{-18,24}},
            textColor={192,192,192},
            textString="tan")}),
      Documentation(info="<html>
<p>
This function returns y = tan(u), with -&infin; &lt; u &lt; &infin;
(if u is a multiple of (2n-1)*pi/2, y = tan(u) is +/- infinity).
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/tan.png\">
</p>
</html>"));
  end tan;

  function asin "Inverse sine (-1 <= u <= 1)"
    extends Modelica.Math.Icons.AxisCenter;
    input Real u "Independent variable";
    output Modelica.Units.SI.Angle y "Dependent variable y=asin(u)";

  external "builtin" y = asin(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,0},{68,0}}, color={192,192,192}),
          Polygon(
            points={{90,0},{68,8},{68,-8},{90,0}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-80},{-79.2,-72.8},{-77.6,-67.5},{-73.6,-59.4},{-66.3,
                -49.8},{-53.5,-37.3},{-30.2,-19.7},{37.4,24.8},{57.5,40.8},{68.7,
                52.7},{75.2,62.2},{77.6,67.5},{80,80}}),
          Text(
            extent={{-88,78},{-16,30}},
            textColor={192,192,192},
            textString="asin")}),
      Documentation(info="<html>
<p>
This function returns y = asin(u), with -1 &le; u &le; +1:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/asin.png\">
</p>
</html>"));
  end asin;

  function cosh "Hyperbolic cosine"
    extends Modelica.Math.Icons.AxisCenter;
    input Real u "Independent variable";
    output Real y "Dependent variable y=cosh(u)";

  external "builtin" y = cosh(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,-86.083},{68,-86.083}}, color={192,192,192}),
          Polygon(
            points={{90,-86.083},{68,-78.083},{68,-94.083},{90,-86.083}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,80},{-77.6,61.1},{-74.4,39.3},{-71.2,20.7},{-67.1,
                1.29},{-63.1,-14.6},{-58.3,-29.8},{-52.7,-43.5},{-46.2,-55.1},{-39,
                -64.3},{-30.2,-71.7},{-18.9,-77.1},{-4.42,-79.9},{10.9,-79.1},{
                23.7,-75.2},{34.2,-68.7},{42.2,-60.6},{48.6,-51.2},{54.3,-40},{
                59.1,-27.5},{63.1,-14.6},{67.1,1.29},{71.2,20.7},{74.4,39.3},{
                77.6,61.1},{80,80}}),
          Text(
            extent={{4,66},{66,20}},
            textColor={192,192,192},
            textString="cosh")}),
      Documentation(info="<html>
<p>
This function returns y = cosh(u), with -&infin; &lt; u &lt; &infin;:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/cosh.png\">
</p>
</html>"));
  end cosh;

  function tanh "Hyperbolic tangent"
    extends Modelica.Math.Icons.AxisCenter;
    input Real u "Independent variable";
    output Real y "Dependent variable y=tanh(u)";

  external "builtin" y = tanh(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,0},{68,0}}, color={192,192,192}),
          Polygon(
            points={{90,0},{68,8},{68,-8},{90,0}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-80},{-47.8,-78.7},{-35.8,-75.7},{-27.7,-70.6},{-22.1,
                -64.2},{-17.3,-55.9},{-12.5,-44.3},{-7.64,-29.2},{-1.21,-4.82},{
                6.83,26.3},{11.7,42},{16.5,54.2},{21.3,63.1},{26.9,69.9},{34.2,75},
                {45.4,78.4},{72,79.9},{80,80}}),
          Text(
            extent={{-88,72},{-16,24}},
            textColor={192,192,192},
            textString="tanh")}),
      Documentation(info="<html>
<p>
This function returns y = tanh(u), with -&infin; &lt; u &lt; &infin;:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/tanh.png\">
</p>
</html>"));
  end tanh;

  function exp "Exponential, base e"
    extends Modelica.Math.Icons.AxisCenter;
    input Real u "Independent variable";
    output Real y "Dependent variable y=exp(u)";

  external "builtin" y = exp(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,-80.3976},{68,-80.3976}}, color={192,192,192}),
          Polygon(
            points={{90,-80.3976},{68,-72.3976},{68,-88.3976},{90,-80.3976}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-80},{-31,-77.9},{-6.03,-74},{10.9,-68.4},{23.7,-61},
                {34.2,-51.6},{43,-40.3},{50.3,-27.8},{56.7,-13.5},{62.3,2.23},{
                67.1,18.6},{72,38.2},{76,57.6},{80,80}}),
          Text(
            extent={{-86,50},{-14,2}},
            textColor={192,192,192},
            textString="exp")}),
      Documentation(info="<html>
<p>
This function returns y = exp(u), with -&infin; &lt; u &lt; &infin;:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/exp.png\">
</p>
</html>"));
  end exp;

  function log "Natural (base e) logarithm (u shall be > 0)"
    extends Modelica.Math.Icons.AxisLeft;
    input Real u "Independent variable";
    output Real y "Dependent variable y=ln(u)";

  external "builtin" y = log(u);
    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}}), graphics={
          Line(points={{-90,0},{68,0}}, color={192,192,192}),
          Polygon(
            points={{90,0},{68,8},{68,-8},{90,0}},
            lineColor={192,192,192},
            fillColor={192,192,192},
            fillPattern=FillPattern.Solid),
          Line(points={{-80,-80},{-79.2,-50.6},{-78.4,-37},{-77.6,-28},{-76.8,-21.3},
                {-75.2,-11.4},{-72.8,-1.31},{-69.5,8.08},{-64.7,17.9},{-57.5,28},
                {-47,38.1},{-31.8,48.1},{-10.1,58},{22.1,68},{68.7,78.1},{80,80}}),
          Text(
            extent={{-6,-24},{66,-72}},
            textColor={192,192,192},
            textString="log")}),
      Documentation(info="<html>
<p>
This function returns y = log(10) (the natural logarithm of u),
with u &gt; 0:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/Math/log.png\">
</p>
</html>"));
  end log;
  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},
            {100,100}}), graphics={Line(points={{-80,0},{-68.7,34.2},{-61.5,53.1},
              {-55.1,66.4},{-49.4,74.6},{-43.8,79.1},{-38.2,79.8},{-32.6,76.6},{
              -26.9,69.7},{-21.3,59.4},{-14.9,44.1},{-6.83,21.2},{10.1,-30.8},{17.3,
              -50.2},{23.7,-64.2},{29.3,-73.1},{35,-78.4},{40.6,-80},{46.2,-77.6},
              {51.9,-71.5},{57.5,-61.9},{63.9,-47.2},{72,-24.8},{80,0}}, color={
              0,0,0}, smooth=Smooth.Bezier)}), Documentation(info="<html>
<p>
This package contains <strong>basic mathematical functions</strong> (such as sin(..)),
as well as functions operating on
<a href=\"modelica://Modelica.Math.Vectors\">vectors</a>,
<a href=\"modelica://Modelica.Math.Matrices\">matrices</a>,
<a href=\"modelica://Modelica.Math.Nonlinear\">nonlinear functions</a>, and
<a href=\"modelica://Modelica.Math.BooleanVectors\">Boolean vectors</a>.
</p>

<h4>Main Authors</h4>
<p><a href=\"http://www.robotic.dlr.de/Martin.Otter/\"><strong>Martin Otter</strong></a>
and <strong>Marcus Baur</strong><br>
Deutsches Zentrum f&uuml;r Luft- und Raumfahrt e.V. (DLR)<br>
Institut f&uuml;r Systemdynamik und Regelungstechnik (DLR-SR)<br>
Forschungszentrum Oberpfaffenhofen<br>
D-82234 Wessling<br>
Germany<br>
email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a>
</p>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>June 22, 2019</em>
       by Thomas Beutlich: Functions tempInterpol1/tempInterpol2 moved to ObsoleteModelica4</li>
<li><em>August 24, 2016</em>
       by Christian Kral: added wrapAngle</li>
<li><em>October 21, 2002</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>
       and Christian Schweiger:<br>
       Function tempInterpol2 added.</li>
<li><em>Oct. 24, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Icons for icon and diagram level introduced.</li>
<li><em>June 30, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized.</li>
</ul>
</html>"));
  end Math;

  package Utilities "Library of utility functions dedicated to scripting (operating on files, streams, strings, system)"
    extends Modelica.Icons.UtilitiesPackage;

    package Streams "Read from files and write to files"
      extends Modelica.Icons.FunctionsPackage;

      pure function error "Print error message and cancel all actions - in case of an unrecoverable error"
        extends Modelica.Icons.Function;
        input String string "String to be printed to error message window";
        external "C" ModelicaError(string) annotation(IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaUtilities.h\"", Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Streams.<strong>error</strong>(string);
</pre></blockquote>
<h4>Description</h4>
<p>
In case of an unrecoverable error (i.e., if the solver is unable to recover from the error),
print the string \"string\" as error message and cancel all actions.
This function is semantically equivalent with the built-in function <strong>assert</strong> if called with the (default) <strong>AssertionLevel.error</strong>.
Line breaks are characterized by \"\\n\" in the string.
</p>
<h4>Example</h4>
<blockquote><pre>
Streams.error(\"x (= \" + String(x) + \")\\nhas to be in the range 0 .. 1\");
</pre></blockquote>
<h4>See also</h4>
<p>
<a href=\"modelica://Modelica.Utilities.Streams\">Streams</a>,
<a href=\"modelica://Modelica.Utilities.Streams.print\">Streams.print</a>,
<a href=\"modelica://ModelicaReference.Operators.'assert()'\">ModelicaReference.Operators.'assert()'</a>
<a href=\"modelica://ModelicaReference.Operators.'String()'\">ModelicaReference.Operators.'String()'</a>
</p>
</html>"));
      end error;
      annotation (
        Documentation(info="<html>
<h4>Library content</h4>
<p>
Package <strong>Streams</strong> contains functions to input and output strings
to a message window or on files, as well as reading matrices from file
and writing matrices to file. Note that a string is interpreted
and displayed as html text (e.g., with print(..) or error(..))
if it is enclosed with the Modelica html quotation, e.g.,
</p>
<blockquote><p>
string = \"&lt;html&gt; first line &lt;br&gt; second line &lt;/html&gt;\".
</p></blockquote>
<p>
It is a quality of implementation, whether (a) all tags of html are supported
or only a subset, (b) how html tags are interpreted if the output device
does not allow to display formatted text.
</p>
<p>
In the table below an example call to every function is given:
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><th><strong><em>Function/type</em></strong></th><th><strong><em>Description</em></strong></th></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.print\">print</a>(string)<br>
          <a href=\"modelica://Modelica.Utilities.Streams.print\">print</a>(string,fileName)</td>
      <td> Print string \"string\" or vector of strings to message window or on
           file \"fileName\".</td>
  </tr>
  <tr><td>stringVector =
         <a href=\"modelica://Modelica.Utilities.Streams.readFile\">readFile</a>(fileName)</td>
      <td> Read complete text file and return it as a vector of strings.</td>
  </tr>
  <tr><td>(string, endOfFile) =
         <a href=\"modelica://Modelica.Utilities.Streams.readLine\">readLine</a>(fileName, lineNumber)</td>
      <td>Returns from the file the content of line lineNumber.</td>
  </tr>
  <tr><td>lines =
         <a href=\"modelica://Modelica.Utilities.Streams.countLines\">countLines</a>(fileName)</td>
      <td>Returns the number of lines in a file.</td>
  </tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.error\">error</a>(string)</td>
      <td> Print error message \"string\" to message window
           and cancel all actions</td>
  </tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.close\">close</a>(fileName)</td>
      <td> Close file if it is still open. Ignore call if
           file is already closed or does not exist. </td>
  </tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.readMatrixSize\">readMatrixSize</a>(fileName, matrixName)</td>
      <td> Read dimensions of a Real matrix from a MATLAB MAT file. </td></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.readRealMatrix\">readRealMatrix</a>(fileName, matrixName, nrow, ncol)</td>
      <td> Read a Real matrix from a MATLAB MAT file. </td></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Streams.writeRealMatrix\">writeRealMatrix</a>(fileName, matrixName, matrix, append, format)</td>
      <td> Write Real matrix to a MATLAB MAT file. </td></tr>
</table>
<p>
Use functions <strong>scanXXX</strong> from package
<a href=\"modelica://Modelica.Utilities.Strings\">Strings</a>
to parse a string.
</p>
<p>
If Real, Integer or Boolean values shall be printed
or used in an error message, they have to be first converted
to strings with the builtin operator
<a href=\"modelica://ModelicaReference.Operators.'String()'\">ModelicaReference.Operators.'String()'</a>(...).
Example:
</p>
<blockquote><pre>
<strong>if</strong> x &lt; 0 <strong>or</strong> x &gt; 1 <strong>then</strong>
   Streams.error(\"x (= \" + String(x) + \") has to be in the range 0 .. 1\");
<strong>end if</strong>;
</pre></blockquote>
</html>"));
    end Streams;

    package Strings "Operations on strings"
      extends Modelica.Icons.FunctionsPackage;

      pure function length "Return length of string"
        extends Modelica.Icons.Function;
        input String string;
        output Integer result "Number of characters of string";
      external "C" result = ModelicaStrings_length(string) annotation(IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStrings.h\"", Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Strings.<strong>length</strong>(string);
</pre></blockquote>
<h4>Description</h4>
<p>
Returns the number of characters of \"string\".
</p>
</html>"));
      end length;

      pure function compare "Compare two strings lexicographically"
        extends Modelica.Icons.Function;
        input String string1;
        input String string2;
        input Boolean caseSensitive=true "= false, if case of letters is ignored";
        output Modelica.Utilities.Types.Compare result "Result of comparison";
      external "C" result = ModelicaStrings_compare(string1, string2, caseSensitive) annotation(IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStrings.h\"", Library="ModelicaExternalC");
        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
result = Strings.<strong>compare</strong>(string1, string2);
result = Strings.<strong>compare</strong>(string1, string2, caseSensitive=true);
</pre></blockquote>
<h4>Description</h4>
<p>
Compares two strings. If the optional argument caseSensitive=false,
upper case letters are treated as if they would be lower case letters.
The result of the comparison is returned as:
</p>
<blockquote><pre>
result = Modelica.Utilities.Types.Compare.Less     // string1 &lt; string2
       = Modelica.Utilities.Types.Compare.Equal    // string1 = string2
       = Modelica.Utilities.Types.Compare.Greater  // string1 &gt; string2
</pre></blockquote>
<p>
Comparison is with regards to lexicographical order,
e.g., \"a\" &lt; \"b\";
</p>
</html>"));
      end compare;

      function isEqual "Determine whether two strings are identical"
        extends Modelica.Icons.Function;
        input String string1;
        input String string2;
        input Boolean caseSensitive=true
          "= false, if lower and upper case are ignored for the comparison";
        output Boolean identical "True, if string1 is identical to string2";
      algorithm
        identical :=compare(string1, string2, caseSensitive) == Types.Compare.Equal;
        annotation (
      Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Strings.<strong>isEqual</strong>(string1, string2);
Strings.<strong>isEqual</strong>(string1, string2, caseSensitive=true);
</pre></blockquote>
<h4>Description</h4>
<p>
Compare whether two strings are identical,
optionally ignoring case.
</p>
</html>"));
      end isEqual;

      function isEmpty
        "Return true if a string is empty (has only white space characters)"
        extends Modelica.Icons.Function;
        input String string;
        output Boolean result "True, if string is empty";
    protected
        Integer nextIndex;
        Integer len;
      algorithm
        nextIndex := Strings.Advanced.skipWhiteSpace(string);
        len := Strings.length(string);
        if len < 1 or nextIndex > len then
          result := true;
        else
          result := false;
        end if;

        annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Strings.<strong>isEmpty</strong>(string);
</pre></blockquote>
<h4>Description</h4>
<p>
Returns true if the string has no characters or if the string consists
only of white space characters. Otherwise, false is returned.
</p>

<h4>Example</h4>
<blockquote><pre>
isEmpty(\"\");       // returns true
isEmpty(\"   \");    // returns true
isEmpty(\"  abc\");  // returns false
isEmpty(\"a\");      // returns false
</pre></blockquote>
</html>"));
      end isEmpty;

      package Advanced "Advanced scanning functions"
        extends Modelica.Icons.FunctionsPackage;

        pure function skipWhiteSpace "Scan white space"
          extends Modelica.Icons.Function;
          input String string;
          input Integer startIndex(min=1)=1;
          output Integer nextIndex;
          external "C" nextIndex = ModelicaStrings_skipWhiteSpace(string, startIndex) annotation(IncludeDirectory="modelica://Modelica/Resources/C-Sources", Include="#include \"ModelicaStrings.h\"", Library="ModelicaExternalC");
          annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
nextIndex = <strong>skipWhiteSpace</strong>(string, startIndex);
</pre></blockquote>
<h4>Description</h4>
<p>
Starts scanning of \"string\" at position \"startIndex\" and
skips white space. The function returns nextIndex = index of character
of the first non white space character.
</p>
<h4>See also</h4>
<a href=\"modelica://Modelica.Utilities.Strings.Advanced\">Strings.Advanced</a>.
</html>"));
        end skipWhiteSpace;
        annotation (Documentation(info="<html>
<h4>Library content</h4>
<p>
Package <strong>Strings.Advanced</strong> contains basic scanning
functions. These functions should be <strong>not called</strong> directly, because
it is much simpler to utilize the higher level functions \"Strings.scanXXX\".
The functions of the \"Strings.Advanced\" library provide
the basic interface in order to implement the higher level
functions in package \"Strings\".
</p>
<p>
Library \"Advanced\" provides the following functions:
</p>
<blockquote><pre>
(nextIndex, realNumber)    = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanReal\">scanReal</a>        (string, startIndex, unsigned=false);
(nextIndex, integerNumber) = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanInteger\">scanInteger</a>     (string, startIndex, unsigned=false);
(nextIndex, string2)       = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanString\">scanString</a>      (string, startIndex);
(nextIndex, identifier)    = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanIdentifier\">scanIdentifier</a>  (string, startIndex);
 nextIndex                 = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.skipWhiteSpace\">skipWhiteSpace</a>  (string, startIndex);
 nextIndex                 = <a href=\"modelica://Modelica.Utilities.Strings.Advanced.skipLineComments\">skipLineComments</a>(string, startIndex);
</pre></blockquote>
<p>
All functions perform the following actions:
</p>
<ol>
<li> Scanning starts at character position \"startIndex\" of
     \"string\" (startIndex has a default of 1).</li>
<li> First, white space is skipped, such as blanks (\" \"), tabs (\"\\t\"), or newline (\"\\n\")</li>
<li> Afterwards, the required token is scanned.</li>
<li> If successful, on return nextIndex = index of character
     directly after the found token and the token value is returned
     as second output argument.<br>
     If not successful, on return nextIndex = startIndex.
     </li>
</ol>
<p>
The following additional rules apply for the scanning:
</p>
<ul>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanReal\">scanReal</a>:<br>
     Scans a full number including one optional leading \"+\" or \"-\" (if unsigned=false)
     according to the Modelica grammar. For example, \"+1.23e-5\", \"0.123\" are
     Real numbers, but \".1\" is not.
     Note, an Integer number, such as \"123\" is also treated as a Real number.<br>&nbsp;</li>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanInteger\">scanInteger</a>:<br>
     Scans an Integer number including one optional leading \"+\"
     or \"-\" (if unsigned=false) according to the Modelica (and C/C++) grammar.
     For example, \"+123\", \"20\" are Integer numbers.
     Note, a Real number, such as \"123.4\" is not an Integer and
     scanInteger returns nextIndex = startIndex.<br>&nbsp;</li>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanString\">scanString</a>:<br>
     Scans a String according to the Modelica (and C/C++) grammar, e.g.,
     \"This is a \"string\"\" is a valid string token.<br>&nbsp;</li>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.scanIdentifier\">scanIdentifier</a>:<br>
     Scans a Modelica identifier, i.e., the identifier starts either
     with a letter, followed by letters, digits or \"_\".
     For example, \"w_rel\", \"T12\".<br>&nbsp;</li>
<li> Function <a href=\"modelica://Modelica.Utilities.Strings.Advanced.skipLineComments\">skipLineComments</a><br>
     Skips white space and Modelica (C/C++) line comments iteratively.
     A line comment starts with \"//\" and ends either with an
     end-of-line (\"\\n\") or the end of the \"string\".</li>
</ul>
</html>"));
      end Advanced;
      annotation (
        Documentation(info="<html>
<h4>Library content</h4>
<p>
Package <strong>Strings</strong> contains functions to manipulate strings.
</p>
<p>
In the table below an example
call to every function is given using the <strong>default</strong> options.
</p>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr><th><strong><em>Function</em></strong></th><th><strong><em>Description</em></strong></th></tr>
  <tr><td>len = <a href=\"modelica://Modelica.Utilities.Strings.length\">length</a>(string)</td>
      <td>Returns length of string</td></tr>
  <tr><td>string2 = <a href=\"modelica://Modelica.Utilities.Strings.substring\">substring</a>(string1,startIndex,endIndex)
       </td>
      <td>Returns a substring defined by start and end index</td></tr>
  <tr><td>result = <a href=\"modelica://Modelica.Utilities.Strings.repeat\">repeat</a>(n)<br>
 result = <a href=\"modelica://Modelica.Utilities.Strings.repeat\">repeat</a>(n,string)</td>
      <td>Repeat a blank or a string n times.</td></tr>
  <tr><td>result = <a href=\"modelica://Modelica.Utilities.Strings.compare\">compare</a>(string1, string2)</td>
      <td>Compares two substrings with regards to alphabetical order</td></tr>
  <tr><td>identical =
<a href=\"modelica://Modelica.Utilities.Strings.isEqual\">isEqual</a>(string1,string2)</td>
      <td>Determine whether two strings are identical</td></tr>
  <tr><td>result = <a href=\"modelica://Modelica.Utilities.Strings.count\">count</a>(string,searchString)</td>
      <td>Count the number of occurrences of a string</td></tr>
  <tr>
<td>index = <a href=\"modelica://Modelica.Utilities.Strings.find\">find</a>(string,searchString)</td>
      <td>Find first occurrence of a string in another string</td></tr>
<tr>
<td>index = <a href=\"modelica://Modelica.Utilities.Strings.findLast\">findLast</a>(string,searchString)</td>
      <td>Find last occurrence of a string in another string</td></tr>
  <tr><td>string2 = <a href=\"modelica://Modelica.Utilities.Strings.replace\">replace</a>(string,searchString,replaceString)</td>
      <td>Replace one or all occurrences of a string</td></tr>
  <tr><td>stringVector2 = <a href=\"modelica://Modelica.Utilities.Strings.sort\">sort</a>(stringVector1)</td>
      <td>Sort vector of strings in alphabetic order</td></tr>
  <tr><td>hash = <a href=\"modelica://Modelica.Utilities.Strings.hashString\">hashString</a>(string)</td>
      <td>Create a hash value of a string</td></tr>
  <tr><td>(token, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanToken\">scanToken</a>(string,startIndex)</td>
      <td>Scan for a token (Real/Integer/Boolean/String/Identifier/Delimiter/NoToken)</td></tr>
  <tr><td>(number, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanReal\">scanReal</a>(string,startIndex)</td>
      <td>Scan for a Real constant</td></tr>
  <tr><td>(number, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanInteger\">scanInteger</a>(string,startIndex)</td>
      <td>Scan for an Integer constant</td></tr>
  <tr><td>(boolean, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanBoolean\">scanBoolean</a>(string,startIndex)</td>
      <td>Scan for a Boolean constant</td></tr>
  <tr><td>(string2, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanString\">scanString</a>(string,startIndex)</td>
      <td>Scan for a String constant</td></tr>
  <tr><td>(identifier, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanIdentifier\">scanIdentifier</a>(string,startIndex)</td>
      <td>Scan for an identifier</td></tr>
  <tr><td>(delimiter, index) = <a href=\"modelica://Modelica.Utilities.Strings.scanDelimiter\">scanDelimiter</a>(string,startIndex)</td>
      <td>Scan for delimiters</td></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Strings.scanNoToken\">scanNoToken</a>(string,startIndex)</td>
      <td>Check that remaining part of string consists solely of<br>
          white space or line comments (\"// ...\\n\").</td></tr>
  <tr><td><a href=\"modelica://Modelica.Utilities.Strings.syntaxError\">syntaxError</a>(string,index,message)</td>
      <td> Print a \"syntax error message\" as well as a string and the<br>
           index at which scanning detected an error</td></tr>
</table>
<p>
The functions \"compare\", \"isEqual\", \"count\", \"find\", \"findLast\", \"replace\", \"sort\"
have the optional
input argument <strong>caseSensitive</strong> with default <strong>true</strong>.
If <strong>false</strong>, the operation is carried out without taking
into account whether a character is upper or lower case.
</p>
</html>"));
    end Strings;

    package Types "Type definitions used in package Modelica.Utilities"
      extends Modelica.Icons.TypesPackage;

      type Compare = enumeration(
          Less "String 1 is lexicographically less than string 2",
          Equal "String 1 is identical to string 2",
          Greater "String 1 is lexicographically greater than string 2")
        "Enumeration defining comparison of two strings";
      annotation (Documentation(info="<html>
<p>
This package contains type definitions used in Modelica.Utilities.
</p>

</html>"));
    end Types;
      annotation (
  Documentation(info="<html>
<p>
This package contains Modelica <strong>functions</strong> that are
especially suited for <strong>scripting</strong>. The functions might
be used to work with strings, read data from file, write data
to file or copy, move and remove files.
</p>
<p>
For an introduction, have especially a look at:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Utilities.UsersGuide\">Modelica.Utilities.User's Guide</a>
     discusses the most important aspects of this library.</li>
<li> <a href=\"modelica://Modelica.Utilities.Examples\">Modelica.Utilities.Examples</a>
     contains examples that demonstrate the usage of this library.</li>
</ul>
<p>
The following main sublibraries are available:
</p>
<ul>
<li> <a href=\"modelica://Modelica.Utilities.Files\">Files</a>
     provides functions to operate on files and directories, e.g.,
     to copy, move, remove files.</li>
<li> <a href=\"modelica://Modelica.Utilities.Streams\">Streams</a>
     provides functions to read from files and write to files.</li>
<li> <a href=\"modelica://Modelica.Utilities.Strings\">Strings</a>
     provides functions to operate on strings. E.g.
     substring, find, replace, sort, scanToken.</li>
<li> <a href=\"modelica://Modelica.Utilities.System\">System</a>
     provides functions to interact with the environment.
     E.g., get or set the working directory or environment
     variables and to send a command to the default shell.</li>
</ul>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>"));
  end Utilities;

  package Constants "Library of mathematical constants and constants of nature (e.g., pi, eps, R, sigma)"
    extends Modelica.Icons.Package;
    import Modelica.Units.SI;
    import Modelica.Units.NonSI;

    final constant Real pi=2*Modelica.Math.asin(1.0);

    final constant Real eps=ModelicaServices.Machine.eps
      "Biggest number such that 1.0 + eps = 1.0";

    final constant Real small=ModelicaServices.Machine.small
      "Smallest number such that small and -small are representable on the machine";

    final constant Real inf=ModelicaServices.Machine.inf
      "Biggest Real number such that inf and -inf are representable on the machine";

    final constant SI.Acceleration g_n=9.80665
      "Standard acceleration of gravity on earth";

    final constant Real k(final unit="J/K") = 1.380649e-23
      "Boltzmann constant";

    final constant Real R(final unit="J/(mol.K)") = k*N_A
      "Molar gas constant";

    final constant Real N_A(final unit="1/mol") = 6.02214076e23
      "Avogadro constant";

    final constant NonSI.Temperature_degC T_zero=-273.15
      "Absolute zero temperature";
    annotation (
      Documentation(info="<html>
<p>
This package provides often needed constants from mathematics, machine
dependent constants and constants from nature. The latter constants
(name, value, description) are from the following source (based on the second source):
</p>
<dl>
<dt>Michael Stock, Richard Davis, Estefan&iacute;a de Mirand&eacute;s and Martin J T Milton:</dt>
<dd><strong>The revision of the SI-the result of three decades of progress in metrology</strong> in Metrologia, Volume 56, Number 2.
<a href= \"https://iopscience.iop.org/article/10.1088/1681-7575/ab0013/pdf\">https://iopscience.iop.org/article/10.1088/1681-7575/ab0013/pdf</a>, 2019.
</dd>
</dl>
<dl>
<dt>D B Newell, F Cabiati, J Fischer, K Fujii, S G Karshenboim, H S Margolis , E de Mirand&eacute;s, P J Mohr, F Nez, K Pachucki, T J Quinn, B N Taylor, M Wang, B M Wood and Z Zhang:</dt>
<dd><strong>The CODATA 2017 values of h, e, k, and NA for the revision of the SI</strong> in Metrologia, Volume 55, Number 1.
<a href= \"https://iopscience.iop.org/article/10.1088/1681-7575/aa950a/pdf\">https://iopscience.iop.org/article/10.1088/1681-7575/aa950a/pdf</a>, 2017.
</dd>
</dl>
<p>BIPM is Bureau International des Poids et Mesures (they publish the SI-standard).</p>
<p>CODATA is the Committee on Data for Science and Technology.</p>

<dl>
<dt><strong>Main Author:</strong></dt>
<dd><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a><br>
    Deutsches Zentrum f&uuml;r Luft und Raumfahrt e. V. (DLR)<br>
    Oberpfaffenhofen<br>
    Postfach 1116<br>
    D-82230 We&szlig;ling<br>
    email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a></dd>
</dl>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>Dec 4, 2019</em>
       by Thomas Beutlich:<br>
       Constant G updated according to 2018 CODATA value.</li>
<li><em>Mar 25, 2019</em>
       by Hans Olsson:<br>
       Constants updated according to 2017 CODATA values and new SI-standard.</li>
<li><em>Nov 4, 2015</em>
       by Thomas Beutlich:<br>
       Constants updated according to 2014 CODATA values.</li>
<li><em>Nov 8, 2004</em>
       by Christian Schweiger:<br>
       Constants updated according to 2002 CODATA values.</li>
<li><em>Dec 9, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Constants updated according to 1998 CODATA values. Using names, values
       and description text from this source. Included magnetic and
       electric constant.</li>
<li><em>Sep 18, 1999</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Constants eps, inf, small introduced.</li>
<li><em>Nov 15, 1997</em>
       by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>
       Realized.</li>
</ul>
</html>"),
      Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
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

  package Icons "Library of icons"
    extends Icons.Package;

    partial package Package "Icon for standard packages"

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Rectangle(
              lineColor={200,200,200},
              fillColor={248,248,248},
              fillPattern=FillPattern.HorizontalCylinder,
              extent={{-100.0,-100.0},{100.0,100.0}},
              radius=25.0),
            Rectangle(
              lineColor={128,128,128},
              extent={{-100.0,-100.0},{100.0,100.0}},
              radius=25.0)}), Documentation(info="<html>
<p>Standard package icon.</p>
</html>"));
    end Package;

    partial package BasesPackage "Icon for packages containing base classes"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Ellipse(
              extent={{-30.0,-30.0},{30.0,30.0}},
              lineColor={128,128,128},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid)}),
                                Documentation(info="<html>
<p>This icon shall be used for a package/library that contains base models and classes, respectively.</p>
</html>"));
    end BasesPackage;

    partial package VariantsPackage "Icon for package containing variants"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                {100,100}}), graphics={
            Ellipse(
              origin={10.0,10.0},
              fillColor={76,76,76},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-80.0,-80.0},{-20.0,-20.0}}),
            Ellipse(
              origin={10.0,10.0},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{0.0,-80.0},{60.0,-20.0}}),
            Ellipse(
              origin={10.0,10.0},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{0.0,0.0},{60.0,60.0}}),
            Ellipse(
              origin={10.0,10.0},
              lineColor={128,128,128},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              extent={{-80.0,0.0},{-20.0,60.0}})}),
                                Documentation(info="<html>
<p>This icon shall be used for a package/library that contains several variants of one component.</p>
</html>"));
    end VariantsPackage;

    partial package InterfacesPackage "Icon for packages containing interfaces"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Polygon(origin={20.0,0.0},
              lineColor={64,64,64},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid,
              points={{-10.0,70.0},{10.0,70.0},{40.0,20.0},{80.0,20.0},{80.0,-20.0},{40.0,-20.0},{10.0,-70.0},{-10.0,-70.0}}),
            Polygon(fillColor={102,102,102},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-100.0,20.0},{-60.0,20.0},{-30.0,70.0},{-10.0,70.0},{-10.0,-70.0},{-30.0,-70.0},{-60.0,-20.0},{-100.0,-20.0}})}),
                                Documentation(info="<html>
<p>This icon indicates packages containing interfaces.</p>
</html>"));
    end InterfacesPackage;

    partial package SourcesPackage "Icon for packages containing sources"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Polygon(origin={23.3333,0.0},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-23.333,30.0},{46.667,0.0},{-23.333,-30.0}}),
            Rectangle(
              fillColor = {128,128,128},
              pattern = LinePattern.None,
              fillPattern = FillPattern.Solid,
              extent = {{-70,-4.5},{0,4.5}})}),
                                Documentation(info="<html>
<p>This icon indicates a package which contains sources.</p>
</html>"));
    end SourcesPackage;

    partial package SensorsPackage "Icon for packages containing sensors"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Ellipse(origin={0.0,-30.0},
              fillColor={255,255,255},
              extent={{-90.0,-90.0},{90.0,90.0}},
              startAngle=20.0,
              endAngle=160.0),
            Ellipse(origin={0.0,-30.0},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-20.0,-20.0},{20.0,20.0}}),
            Line(origin={0.0,-30.0},
              points={{0.0,60.0},{0.0,90.0}}),
            Ellipse(origin={-0.0,-30.0},
              fillColor={64,64,64},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-10.0,-10.0},{10.0,10.0}}),
            Polygon(
              origin={-0.0,-30.0},
              rotation=-35.0,
              fillColor={64,64,64},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-7.0,0.0},{-3.0,85.0},{0.0,90.0},{3.0,85.0},{7.0,0.0}})}),
                                Documentation(info="<html>
<p>This icon indicates a package containing sensors.</p>
</html>"));
    end SensorsPackage;

    partial package UtilitiesPackage "Icon for utility packages"
      extends Modelica.Icons.Package;
       annotation (Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
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
        points={{-15.0,87.273},{15.0,87.273},{20.0,82.273},{20.0,27.273},{10.0,17.273},{10.0,7.273},{20.0,2.273},{20.0,-2.727},{5.0,-2.727},{5.0,-77.727},{10.0,-87.727},{5.0,-112.727},{-5.0,-112.727},{-10.0,-87.727},{-5.0,-77.727},{-5.0,-2.727},{-20.0,-2.727},{-20.0,2.273},{-10.0,7.273},{-10.0,17.273},{-20.0,27.273},{-20.0,82.273}})}),
      Documentation(info="<html>
<p>This icon indicates a package containing utility classes.</p>
</html>"));
    end UtilitiesPackage;

    partial package TypesPackage "Icon for packages containing type definitions"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-12.167,-23},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{12.167,65},{14.167,93},{36.167,89},{24.167,20},{4.167,-30},
                  {14.167,-30},{24.167,-30},{24.167,-40},{-5.833,-50},{-15.833,
                  -30},{4.167,20},{12.167,65}},
              smooth=Smooth.Bezier), Polygon(
              origin={2.7403,1.6673},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{49.2597,22.3327},{31.2597,24.3327},{7.2597,18.3327},{-26.7403,
                10.3327},{-46.7403,14.3327},{-48.7403,6.3327},{-32.7403,0.3327},{-6.7403,
                4.3327},{33.2597,14.3327},{49.2597,14.3327},{49.2597,22.3327}},
              smooth=Smooth.Bezier)}));
    end TypesPackage;

    partial package FunctionsPackage "Icon for packages containing functions"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
              Text(
                textColor={128,128,128},
                extent={{-90,-90},{90,90}},
                textString="f")}));
    end FunctionsPackage;

    partial package IconsPackage "Icon for packages containing icons"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-8.167,-17},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-15.833,20.0},{-15.833,30.0},{14.167,40.0},{24.167,20.0},{
                  4.167,-30.0},{14.167,-30.0},{24.167,-30.0},{24.167,-40.0},{-5.833,
                  -50.0},{-15.833,-30.0},{4.167,20.0},{-5.833,20.0}},
              smooth=Smooth.Bezier), Ellipse(
              origin={-0.5,56.5},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-12.5,-12.5},{12.5,12.5}})}));
    end IconsPackage;

    partial package InternalPackage "Icon for an internal package (indicating that the package should not be directly utilized by user)"
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}), graphics={
          Rectangle(
            lineColor={215,215,215},
            fillColor={255,255,255},
            fillPattern=FillPattern.HorizontalCylinder,
            extent={{-100,-100},{100,100}},
            radius=25),
          Rectangle(
            lineColor={215,215,215},
            extent={{-100,-100},{100,100}},
            radius=25),
          Ellipse(
            extent={{-80,80},{80,-80}},
            lineColor={215,215,215},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid),
          Ellipse(
            extent={{-55,55},{55,-55}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent={{-60,14},{60,-14}},
            lineColor={215,215,215},
            fillColor={215,215,215},
            fillPattern=FillPattern.Solid,
            rotation=45)}),
      Documentation(info="<html>

<p>
This icon shall be used for a package that contains internal classes not to be
directly utilized by a user.
</p>
</html>"));
    end InternalPackage;

    partial package MaterialPropertiesPackage "Icon for package containing property classes"
      extends Modelica.Icons.Package;
      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={
            Ellipse(
              lineColor={102,102,102},
              fillColor={204,204,204},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Sphere,
              extent={{-60.0,-60.0},{60.0,60.0}})}),
                                Documentation(info="<html>
<p>This icon indicates a package that contains properties</p>
</html>"));
    end MaterialPropertiesPackage;

    partial class RoundSensor "Icon representing a round measurement device"

      annotation (
        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Ellipse(
              fillColor={245,245,245},
              fillPattern=FillPattern.Solid,
              extent={{-70.0,-70.0},{70.0,70.0}}),
            Line(points={{0.0,70.0},{0.0,40.0}}),
            Line(points={{22.9,32.8},{40.2,57.3}}),
            Line(points={{-22.9,32.8},{-40.2,57.3}}),
            Line(points={{37.6,13.7},{65.8,23.9}}),
            Line(points={{-37.6,13.7},{-65.8,23.9}}),
            Ellipse(
              lineColor={64,64,64},
              fillColor={255,255,255},
              extent={{-12.0,-12.0},{12.0,12.0}}),
            Polygon(
              rotation=-17.5,
              fillColor={64,64,64},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-5.0,0.0},{-2.0,60.0},{0.0,65.0},{2.0,60.0},{5.0,0.0}}),
            Ellipse(
              fillColor={64,64,64},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-7.0,-7.0},{7.0,7.0}})}),
        Documentation(info="<html>
<p>
This icon is designed for a <strong>rotational sensor</strong> model.
</p>
</html>"));
    end RoundSensor;

    partial function Function "Icon for functions"

      annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
            Text(
              textColor={0,0,255},
              extent={{-150,105},{150,145}},
              textString="%name"),
            Ellipse(
              lineColor = {108,88,49},
              fillColor = {255,215,136},
              fillPattern = FillPattern.Solid,
              extent = {{-100,-100},{100,100}}),
            Text(
              textColor={108,88,49},
              extent={{-90.0,-90.0},{90.0,90.0}},
              textString="f")}),
    Documentation(info="<html>
<p>This icon indicates Modelica functions.</p>
</html>"));
    end Function;

    partial record Record "Icon for records"

      annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), graphics={
            Text(
              textColor={0,0,255},
              extent={{-150,60},{150,100}},
              textString="%name"),
            Rectangle(
              origin={0.0,-25.0},
              lineColor={64,64,64},
              fillColor={255,215,136},
              fillPattern=FillPattern.Solid,
              extent={{-100.0,-75.0},{100.0,75.0}},
              radius=25.0),
            Line(
              points={{-100.0,0.0},{100.0,0.0}},
              color={64,64,64}),
            Line(
              origin={0.0,-50.0},
              points={{-100.0,0.0},{100.0,0.0}},
              color={64,64,64}),
            Line(
              origin={0.0,-25.0},
              points={{0.0,75.0},{0.0,-75.0}},
              color={64,64,64})}), Documentation(info="<html>
<p>
This icon is indicates a record.
</p>
</html>"));
    end Record;

    type TypeReal "Icon for Real types"
        extends Real;
        annotation(Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
              Rectangle(
                lineColor={160,160,164},
                fillColor={160,160,164},
                fillPattern=FillPattern.Solid,
                extent={{-100.0,-100.0},{100.0,100.0}},
                radius=25.0),
              Text(
                textColor={255,255,255},
                extent={{-90.0,-50.0},{90.0,50.0}},
                textString="R")}),Documentation(info="<html>
<p>
This icon is designed for a <strong>Real</strong> type.
</p>
</html>"));
    end TypeReal;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
                -100},{100,100}}), graphics={Polygon(
              origin={-8.167,-17},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              points={{-15.833,20.0},{-15.833,30.0},{14.167,40.0},{24.167,20.0},{
                  4.167,-30.0},{14.167,-30.0},{24.167,-30.0},{24.167,-40.0},{-5.833,
                  -50.0},{-15.833,-30.0},{4.167,20.0},{-5.833,20.0}},
              smooth=Smooth.Bezier), Ellipse(
              origin={-0.5,56.5},
              fillColor={128,128,128},
              pattern=LinePattern.None,
              fillPattern=FillPattern.Solid,
              extent={{-12.5,-12.5},{12.5,12.5}})}), Documentation(info="<html>
<p>This package contains definitions for the graphical layout of components which may be used in different libraries. The icons can be utilized by inheriting them in the desired class using &quot;extends&quot; or by directly copying the &quot;icon&quot; layer.</p>

<h4>Main Authors</h4>

<dl>
<dt><a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a></dt>
    <dd>Deutsches Zentrum fuer Luft und Raumfahrt e.V. (DLR)</dd>
    <dd>Oberpfaffenhofen</dd>
    <dd>Postfach 1116</dd>
    <dd>D-82230 Wessling</dd>
    <dd>email: <a href=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</a></dd>
<dt>Christian Kral</dt>

    <dd>  <a href=\"https://christiankral.net/\">Electric Machines, Drives and Systems</a><br>
</dd>
    <dd>1060 Vienna, Austria</dd>
    <dd>email: <a href=\"mailto:dr.christian.kral@gmail.com\">dr.christian.kral@gmail.com</a></dd>
<dt>Johan Andreasson</dt>
    <dd><a href=\"https://www.modelon.com/\">Modelon AB</a></dd>
    <dd>Ideon Science Park</dd>
    <dd>22370 Lund, Sweden</dd>
    <dd>email: <a href=\"mailto:johan.andreasson@modelon.se\">johan.andreasson@modelon.se</a></dd>
</dl>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>"));
  end Icons;

  package Units "Library of type and unit definitions"
    extends Modelica.Icons.Package;

    package SI "Library of SI unit definitions"
      extends Modelica.Icons.Package;

      type Angle = Real (
          final quantity="Angle",
          final unit="rad",
          displayUnit="deg");

      type Length = Real (final quantity="Length", final unit="m");

      type Height = Length(min=0);

      type Diameter = Length(min=0);

      type Area = Real (final quantity="Area", final unit="m2");

      type Volume = Real (final quantity="Volume", final unit="m3");

      type Time = Real (final quantity="Time", final unit="s");

      type Velocity = Real (final quantity="Velocity", final unit="m/s");

      type Acceleration = Real (final quantity="Acceleration", final unit="m/s2");

      type Mass = Real (
          quantity="Mass",
          final unit="kg",
          min=0);

      type Density = Real (
          final quantity="Density",
          final unit="kg/m3",
          displayUnit="g/cm3",
          min=0.0);

      type Pressure = Real (
          final quantity="Pressure",
          final unit="Pa",
          displayUnit="bar");

      type AbsolutePressure = Pressure (min=0.0, nominal = 1e5);

      type PressureDifference = Pressure;

      type DynamicViscosity = Real (
          final quantity="DynamicViscosity",
          final unit="Pa.s",
          min=0);

      type Energy = Real (final quantity="Energy", final unit="J");

      type Power = Real (final quantity="Power", final unit="W");

      type EnergyFlowRate = Power;

      type EnthalpyFlowRate = Real (final quantity="EnthalpyFlowRate", final unit=
              "W");

      type Efficiency = Real (
          final quantity="Efficiency",
          final unit="1",
          min=0);

      type MassFlowRate = Real (quantity="MassFlowRate", final unit="kg/s");

      type VolumeFlowRate = Real (final quantity="VolumeFlowRate", final unit=
              "m3/s");

      type ThermodynamicTemperature = Real (
          final quantity="ThermodynamicTemperature",
          final unit="K",
          min = 0.0,
          start = 288.15,
          nominal = 300,
          displayUnit="degC")
        "Absolute temperature (use type TemperatureDifference for relative temperatures)" annotation(absoluteValue=true);

      type Temperature = ThermodynamicTemperature;

      type TemperatureDifference = Real (
          final quantity="ThermodynamicTemperature",
          final unit="K") annotation(absoluteValue=false);

      type TemperatureSlope = Real (final quantity="TemperatureSlope",
          final unit="K/s");

      type LinearTemperatureCoefficient = Real(final quantity = "LinearTemperatureCoefficient", final unit="1/K");

      type Compressibility = Real (final quantity="Compressibility", final unit=
              "1/Pa");

      type IsothermalCompressibility = Compressibility;

      type HeatFlowRate = Real (final quantity="Power", final unit="W");

      type ThermalConductivity = Real (final quantity="ThermalConductivity", final unit=
                 "W/(m.K)");

      type CoefficientOfHeatTransfer = Real (final quantity=
              "CoefficientOfHeatTransfer", final unit="W/(m2.K)");

      type ThermalResistance = Real (final quantity="ThermalResistance", final unit=
             "K/W");

      type ThermalConductance = Real (final quantity="ThermalConductance", final unit=
                 "W/K");

      type HeatCapacity = Real (final quantity="HeatCapacity", final unit="J/K");

      type SpecificHeatCapacity = Real (final quantity="SpecificHeatCapacity",
            final unit="J/(kg.K)");

      type RatioOfSpecificHeatCapacities = Real (final quantity=
              "RatioOfSpecificHeatCapacities", final unit="1");

      type SpecificEntropy = Real (final quantity="SpecificEntropy",
                                   final unit="J/(kg.K)");

      type SpecificEnergy = Real (final quantity="SpecificEnergy",
                                  final unit="J/kg");

      type SpecificInternalEnergy = SpecificEnergy;

      type SpecificEnthalpy = SpecificEnergy;

      type DerDensityByEnthalpy = Real (final unit="kg.s2/m5");

      type DerDensityByPressure = Real (final unit="s2/m2");

      type DerDensityByTemperature = Real (final unit="kg/(m3.K)");

      type MolarMass = Real (final quantity="MolarMass", final unit="kg/mol", min=0);

      type MolarVolume = Real (final quantity="MolarVolume", final unit="m3/mol", min=0);

      type MassFraction = Real (final quantity="MassFraction", final unit="1",
                                min=0, max=1);

      type MoleFraction = Real (final quantity="MoleFraction", final unit="1",
                                min = 0, max = 1);

      type ReynoldsNumber = Real (final quantity="ReynoldsNumber", final unit="1");

      type PrandtlNumber = Real (final quantity="PrandtlNumber", final unit="1");
      annotation (Icon(graphics={Text(
              extent={{-80,80},{80,-78}},
              textColor={128,128,128},
              fillColor={128,128,128},
              fillPattern=FillPattern.None,
              fontName="serif",
              textString="SI",
              textStyle={TextStyle.Italic})}),
                                       Documentation(info="<html>
<p>This package provides predefined types based on the international standard
on units.
</p>
<p>
For an introduction to the conventions used in this package, have a look at:
<a href=\"modelica://Modelica.Units.UsersGuide.Conventions\">Conventions</a>.
</p>
</html>"));
    end SI;

    package NonSI "Type definitions of non SI and other units"
      extends Modelica.Icons.Package;

      type Temperature_degC = Real (final quantity="ThermodynamicTemperature",
            final unit="degC")
        "Absolute temperature in degree Celsius (for relative temperature use Modelica.Units.SI.TemperatureDifference)" annotation(absoluteValue=true);

      type Pressure_bar = Real (final quantity="Pressure", final unit="bar")
        "Absolute pressure in bar";
      annotation (Documentation(info="<html>
<p>
This package provides predefined types, such as <strong>Angle_deg</strong> (angle in
degree), <strong>AngularVelocity_rpm</strong> (angular velocity in revolutions per
minute) or <strong>Temperature_degF</strong> (temperature in degree Fahrenheit),
which are in common use but are not part of the international standard on
units according to ISO 31-1992 \"General principles concerning quantities,
units and symbols\" and ISO 1000-1992 \"SI units and recommendations for
the use of their multiples and of certain other units\".</p>
<p>If possible, the types in this package should not be used. Use instead
types of package <code>Modelica.Units.SI</code>. For more information on units, see also
the book of Francois Cardarelli <strong>Scientific Unit Conversion - A
Practical Guide to Metrication</strong> (Springer 1997).</p>
</html>"), Icon(coordinateSystem(extent={{-100,-100},{100,100}}), graphics={Ellipse(
              extent={{-10,10},{10,-10}},
              lineColor={128,128,128},
              fillColor={128,128,128},
              fillPattern=FillPattern.Solid), Ellipse(
              extent={{-60,10},{-40,-10}},
              lineColor={128,128,128},
              fillColor={128,128,128},
              fillPattern=FillPattern.Solid), Ellipse(
              extent={{40,10},{60,-10}},
              lineColor={128,128,128},
              fillColor={128,128,128},
              fillPattern=FillPattern.Solid)}));
    end NonSI;

    package Conversions "Conversion functions to/from non SI units and type definitions of non SI units"
      extends Modelica.Icons.Package;

      function to_degC "Convert from kelvin to degree Celsius"
        extends Modelica.Units.Icons.Conversion;
        input SI.Temperature Kelvin "Value in kelvin";
        output Modelica.Units.NonSI.Temperature_degC Celsius "Value in degree Celsius";
      algorithm
        Celsius := Kelvin + Modelica.Constants.T_zero;
        annotation (Inline=true,Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-20,100},{-100,20}},
                textString="K"), Text(
                extent={{100,-20},{20,-100}},
                textString="degC")}));
      end to_degC;

      function from_degC "Convert from degree Celsius to kelvin"
        extends Modelica.Units.Icons.Conversion;
        input Modelica.Units.NonSI.Temperature_degC Celsius "Value in degree Celsius";
        output SI.Temperature Kelvin "Value in kelvin";
      algorithm
        Kelvin := Celsius - Modelica.Constants.T_zero;
        annotation (Inline=true,Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-20,100},{-100,20}},
                textString="degC"), Text(
                extent={{100,-20},{20,-100}},
                textString="K")}));
      end from_degC;

      function to_bar "Convert from Pascal to bar"
        extends Modelica.Units.Icons.Conversion;
        input SI.Pressure Pa "Value in Pascal";
        output Modelica.Units.NonSI.Pressure_bar bar "Value in bar";
      algorithm
        bar := Pa/1e5;
        annotation (Inline=true,Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={Text(
                extent={{-12,100},{-100,56}},
                textString="Pa"), Text(
                extent={{98,-52},{-4,-100}},
                textString="bar")}));
      end to_bar;
      annotation (Documentation(info="<html>
<p>This package provides conversion functions from the non SI Units
defined in package <code>Modelica.Units.NonSI</code> to the
corresponding SI Units defined in package <code>Modelica.Units.SI</code> and vice
versa. It is recommended to use these functions in the following
way (note, that all functions have one Real input and one Real output
argument):</p>
<blockquote><pre>
<strong>import</strong> Modelica.Units.SI;
<strong>import</strong> Modelica.Units.Conversions.{from_degC, from_deg, from_rpm};
   ...
<strong>parameter</strong> SI.Temperature     T   = from_degC(25);   // convert 25 degree Celsius to kelvin
<strong>parameter</strong> SI.Angle           phi = from_deg(180);   // convert 180 degree to radian
<strong>parameter</strong> SI.AngularVelocity w   = from_rpm(3600);  // convert 3600 revolutions per minutes
                                                   // to radian per seconds
</pre></blockquote>

</html>"),     Icon(graphics={
            Polygon(
              points={{80,0},{20,20},{20,-20},{80,0}},
              lineColor={191,0,0},
              fillColor={191,0,0},
              fillPattern=FillPattern.Solid),
            Line(points={{-80,0},{20,0}}, color={191,0,0})}));
    end Conversions;

    package Icons "Icons for Units"
      extends Modelica.Icons.IconsPackage;

      partial function Conversion "Base icon for conversion functions"

        annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
                  -100},{100,100}}), graphics={
              Rectangle(
                extent={{-100,100},{100,-100}},
                lineColor={191,0,0},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid),
              Line(points={{-90,0},{30,0}}, color={191,0,0}),
              Polygon(
                points={{90,0},{30,20},{30,-20},{90,0}},
                lineColor={191,0,0},
                fillColor={191,0,0},
                fillPattern=FillPattern.Solid),
              Text(
                extent={{-115,155},{115,105}},
                textString="%name",
                textColor={0,0,255})}));
      end Conversion;
    end Icons;
    annotation (Icon(graphics={
        Polygon(
          fillColor = {128,128,128},
          pattern = LinePattern.None,
          fillPattern = FillPattern.Solid,
          points = {{-80,-40},{-80,-40},{-55,50},{-52.5,62.5},{-65,60},{-65,65},{-35,77.5},{-32.5,60},{-50,0},{-50,0},{-30,15},{-20,27.5},{-32.5,27.5},{-32.5,27.5},{-32.5,32.5},{-32.5,32.5},{2.5,32.5},{2.5,32.5},{2.5,27.5},{2.5,27.5},{-7.5,27.5},{-30,7.5},{-30,7.5},{-25,-25},{-17.5,-28.75},{-10,-25},{-5,-26.25},{-5,-32.5},{-16.25,-41.25},{-31.25,-43.75},{-40,-33.75},{-45,-5},{-45,-5},{-52.5,-10},{-52.5,-10},{-60,-40},{-60,-40}},
          smooth = Smooth.Bezier),
        Polygon(
          fillColor = {128,128,128},
          pattern = LinePattern.None,
          fillPattern = FillPattern.Solid,
          points = {{87.5,30},{62.5,30},{62.5,30},{55,33.75},{36.25,35},{16.25,25},{7.5,6.25},{11.25,-7.5},{22.5,-12.5},{22.5,-12.5},{6.25,-22.5},{6.25,-35},{16.25,-38.75},{16.25,-38.75},{21.25,-41.25},{21.25,-41.25},{45,-48.75},{47.5,-61.25},{32.5,-70},{12.5,-65},{7.5,-51.25},{21.25,-41.25},{21.25,-41.25},{16.25,-38.75},{16.25,-38.75},{6.25,-41.25},{-6.25,-50},{-3.75,-68.75},{30,-76.25},{65,-62.5},{63.75,-35},{27.5,-26.25},{22.5,-20},{27.5,-15},{27.5,-15},{30,-7.5},{30,-7.5},{27.5,-2.5},{28.75,11.25},{36.25,27.5},{47.5,30},{53.75,22.5},{51.25,8.75},{45,-6.25},{35,-11.25},{30,-7.5},{30,-7.5},{27.5,-15},{27.5,-15},{43.75,-16.25},{65,-6.25},{72.5,10},{70,20},{70,20},{80,20}},
          smooth = Smooth.Bezier)}), Documentation(info="<html>
<p>This package provides predefined types, such as <em>Mass</em>,
<em>Angle</em>, <em>Time</em>, based on the international standard
on units, e.g.,
</p>

<blockquote><pre>
<strong>type</strong> Angle = Real(<strong>final</strong> quantity = \"Angle\",
                  <strong>final</strong> unit     = \"rad\",
                  displayUnit   = \"deg\");
</pre></blockquote>

<p>
Some of the types are derived SI units that are utilized in package Modelica
(such as ComplexCurrent, which is a complex number where both the real and imaginary
part have the SI unit Ampere).
</p>

<p>
Furthermore, conversion functions from non SI-units to SI-units and vice versa
are provided in subpackage
<a href=\"modelica://Modelica.Units.Conversions\">Conversions</a>.
</p>

<p>
For an introduction how units are used in the Modelica Standard Library
with package Units, have a look at:
<a href=\"modelica://Modelica.Units.UsersGuide.HowToUseUnits\">How to use Units</a>.
</p>

<p>
Copyright &copy; 1998-2020, Modelica Association and contributors
</p>
</html>",   revisions="<html>
<ul>
<li><em>May 25, 2011</em> by Stefan Wischhusen:<br>Added molar units for energy and enthalpy.</li>
<li><em>Jan. 27, 2010</em> by Christian Kral:<br>Added complex units.</li>
<li><em>Dec. 14, 2005</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Add User&#39;s Guide and removed &quot;min&quot; values for Resistance and Conductance.</li>
<li><em>October 21, 2002</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Christian Schweiger:<br>Added new package <strong>Conversions</strong>. Corrected typo <em>Wavelenght</em>.</li>
<li><em>June 6, 2000</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Introduced the following new types<br>type Temperature = ThermodynamicTemperature;<br>types DerDensityByEnthalpy, DerDensityByPressure, DerDensityByTemperature, DerEnthalpyByPressure, DerEnergyByDensity, DerEnergyByPressure<br>Attribute &quot;final&quot; removed from min and max values in order that these values can still be changed to narrow the allowed range of values.<br>Quantity=&quot;Stress&quot; removed from type &quot;Stress&quot;, in order that a type &quot;Stress&quot; can be connected to a type &quot;Pressure&quot;.</li>
<li><em>Oct. 27, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>New types due to electrical library: Transconductance, InversePotential, Damping.</li>
<li><em>Sept. 18, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Renamed from SIunit to SIunits. Subpackages expanded, i.e., the SIunits package, does no longer contain subpackages.</li>
<li><em>Aug 12, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Type &quot;Pressure&quot; renamed to &quot;AbsolutePressure&quot; and introduced a new type &quot;Pressure&quot; which does not contain a minimum of zero in order to allow convenient handling of relative pressure. Redefined BulkModulus as an alias to AbsolutePressure instead of Stress, since needed in hydraulics.</li>
<li><em>June 29, 1999</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a>:<br>Bug-fix: Double definition of &quot;Compressibility&quot; removed and appropriate &quot;extends Heat&quot; clause introduced in package SolidStatePhysics to incorporate ThermodynamicTemperature.</li>
<li><em>April 8, 1998</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Astrid Jaschinski:<br>Complete ISO 31 chapters realized.</li>
<li><em>Nov. 15, 1997</em> by <a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a> and Hubertus Tummescheit:<br>Some chapters realized.</li>
</ul>
</html>"));
  end Units;
annotation (
preferredView="info",
version="4.0.0",
versionDate="2020-06-04",
dateModified = "2020-06-04 11:00:00Z",
revisionId="6626538a2 2020-06-04 19:56:34 +0200",
uses(Complex(version="4.0.0"), ModelicaServices(version="4.0.0")),
conversion(
 from(version={"3.0", "3.0.1", "3.1", "3.2", "3.2.1", "3.2.2", "3.2.3"}, script="modelica://Modelica/Resources/Scripts/Conversion/ConvertModelica_from_3.2.3_to_4.0.0.mos")),
Icon(coordinateSystem(extent={{-100.0,-100.0},{100.0,100.0}}), graphics={
  Polygon(
    origin={-6.9888,20.048},
    pattern=LinePattern.None,
    fillPattern=FillPattern.Solid,
    points={{-93.0112,10.3188},{-93.0112,10.3188},{-73.011,24.6},{-63.011,31.221},{-51.219,36.777},{-39.842,38.629},{-31.376,36.248},{-25.819,29.369},{-24.232,22.49},{-23.703,17.463},{-15.501,25.135},{-6.24,32.015},{3.02,36.777},{15.191,39.423},{27.097,37.306},{32.653,29.633},{35.035,20.108},{43.501,28.046},{54.085,35.19},{65.991,39.952},{77.897,39.688},{87.422,33.338},{91.126,21.696},{90.068,9.525},{86.099,-1.058},{79.749,-10.054},{71.283,-21.431},{62.816,-33.337},{60.964,-32.808},{70.489,-16.14},{77.368,-2.381},{81.072,10.054},{79.749,19.05},{72.605,24.342},{61.758,23.019},{49.587,14.817},{39.003,4.763},{29.214,-6.085},{21.012,-16.669},{13.339,-26.458},{5.401,-36.777},{-1.213,-46.037},{-6.24,-53.446},{-8.092,-52.387},{-0.684,-40.746},{5.401,-30.692},{12.81,-17.198},{19.424,-3.969},{23.658,7.938},{22.335,18.785},{16.514,23.283},{8.047,23.019},{-1.478,19.05},{-11.267,11.113},{-19.734,2.381},{-29.259,-8.202},{-38.519,-19.579},{-48.044,-31.221},{-56.511,-43.392},{-64.449,-55.298},{-72.386,-66.939},{-77.678,-74.612},{-79.53,-74.083},{-71.857,-61.383},{-62.861,-46.037},{-52.278,-28.046},{-44.869,-15.346},{-38.784,-2.117},{-35.344,8.731},{-36.403,19.844},{-42.488,23.813},{-52.013,22.49},{-60.744,16.933},{-68.947,10.054},{-76.884,2.646},{-93.0112,-12.1707},{-93.0112,-12.1707}},
    smooth=Smooth.Bezier),
  Ellipse(
    origin={40.8208,-37.7602},
    fillColor={161,0,4},
    pattern=LinePattern.None,
    fillPattern=FillPattern.Solid,
    extent={{-17.8562,-17.8563},{17.8563,17.8562}})}),
Documentation(info="<html>
<p>
<img src=\"modelica://Modelica/Resources/Images/Logos/Modelica_Libraries.svg\" width=\"250\">
</p>

<p>
The package <strong>Modelica&reg;</strong> is a <strong>standardized</strong> and <strong>free</strong> package
that is developed by the \"<strong>Modelica Association Project - Libraries</strong>\".</p>
<p>
Its development is coordinated with the Modelica&reg; language from the
Modelica Association, see <a href=\"https://www.Modelica.org\">https://www.Modelica.org</a>.
It is also called <strong>Modelica Standard Library</strong>.
It provides model components in many domains that are based on
standardized interface definitions. Some typical examples are shown
in the next figure:
</p>

<p>
<img src=\"modelica://Modelica/Resources/Images/UsersGuide/ModelicaLibraries.png\">
</p>

<p>
For an introduction, have especially a look at:
</p>
<ul>
<li> <a href=\"modelica://Modelica.UsersGuide.Overview\">Overview</a>
  provides an overview of the Modelica Standard Library
  inside the <a href=\"modelica://Modelica.UsersGuide\">User's Guide</a>.</li>
<li><a href=\"modelica://Modelica.UsersGuide.ReleaseNotes\">Release Notes</a>
 summarizes the changes of new versions of this package.</li>
<li> <a href=\"modelica://Modelica.UsersGuide.Contact\">Contact</a>
  lists the contributors of the Modelica Standard Library.</li>
<li> The <strong>Examples</strong> packages in the various libraries, demonstrate
  how to use the components of the corresponding sublibrary.</li>
</ul>

<p>
This version of the Modelica Standard Library consists of
</p>
<ul>
<li><strong>1417</strong> component models and blocks,</li>
<li><strong>512</strong> example models, and</li>
<li><strong>1219</strong> functions</li>
</ul>
<p>
that are directly usable (= number of public, non-partial, non-internal and non-obsolete classes). It is fully compliant
to <a href=\"https://modelica.org/documents/ModelicaSpec34.pdf\">Modelica Specification version 3.4</a>
and it has been tested with Modelica tools from different vendors.
</p>

<p>
<strong>Licensed by the Modelica Association under the 3-Clause BSD License</strong><br>
Copyright &copy; 1998-2020, Modelica Association and <a href=\"modelica://Modelica.UsersGuide.Contact\">contributors</a>.
</p>

<p>
<em>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the 3-Clause BSD license. For license conditions (including the disclaimer of warranty) visit <a href=\"https://modelica.org/licenses/modelica-3-clause-bsd\">https://modelica.org/licenses/modelica-3-clause-bsd</a>.</em>
</p>

<p>
<strong>Modelica&reg;</strong> is a registered trademark of the Modelica Association.
</p>
</html>"));
end Modelica;
