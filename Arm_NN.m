clear
% Node and timing information
N = 25;                 % # of nodal points
duration = 1.0;         % time in sec to complete task
h = duration/(N-1);     % time interval between nodes
dc_time = h*(0:N-1)';   % list of time points (temporal grid)

% Import the OpenSim modeling classes
import org.opensim.modeling.*

% Read in the osim model
osimModel = Model('E:\Arm26\Arm26.osim'); % 

% Initialize the model (this builds the system and initialize the state)
osimState = osimModel.initSystem();

% Get the number of states, coordinates, muscles and controls from the model;
% in this case the number of controls equals the number of muscles

Nstates       = osimModel.getNumStateVariables();
Ncontrols     = osimModel.getNumControls();
Ncoord        = osimModel.getNumCoordinates(); 
model_muscles = osimModel.getMuscles();
Nmuscles      = model_muscles.getSize();

% Auxiliary data 

auxdata.model      = osimModel;
auxdata.time       = dc_time;
auxdata.N          = N;
auxdata.h          = h;
auxdata.Nstates    = Nstates;
auxdata.Ncontrols  = Ncontrols; 
auxdata.Nmuscles   = Nmuscles;
auxdata.Ncoord     = Ncoord;

% Get the names of the states from the model

states_all = cell(Nstates,1);
for i = 1:Nstates
   states_all(i,1) = cell(osimModel.getStateVariableNames().getitem(i-1));
end

% Get the names of the controls/muscles from the model (same in this case)

Muscles = osimModel.getMuscles();  
controls_all = cell(Ncontrols,1);
for i = 1:Ncontrols
   currentMuscle = Muscles.get(i-1);
   controls_all(i,1) = cell(currentMuscle.getName());
end



[file_input, pathname] = uigetfile({'*.sto', 'OpenSim States Files (*.sto)'}, ...
                         'Select the initial states file','MultiSelect', 'off');
temp_s = importdata(strcat(pathname,file_input)); % import states data
old_time_s = temp_s.data(:,1);         % time 
old_data_s = temp_s.data(:,2:end);     % the 6 states
old_text_s = temp_s.textdata(7,2:end); % states names

% Arrange the initial guess by nodes and states

x0_temp = zeros(N,Nstates);  % pre-allocate space
chk_counter =0;
for j = 1:size(states_all,1) 
    for k = 1:size(old_data_s,2)
        if strcmp(old_text_s(k),states_all(j)) == 1,
            % interpolate initial guess to the defined temporal grid
            x0_temp(:,j) = interp1(old_time_s,old_data_s(:,k),dc_time);
            chk_counter = chk_counter+1;
        end
    end
end

if (chk_counter ~= size(states_all,1))
    disp('Inconsistent number of states for inital guess!');
    disp('Check the input file...');
    return
end


% Arrange the initial guess into a column vector
% [ [pos(t_0) ... pos(t_N)], [vel(t_0) ... vel(t_N)], ... etc]

for i = 1:Nstates,
    X0(N*(i-1)+1:N*i,:) = x0_temp(:,i);
end

% Load the file that contain the initial guess for the controls (excitations)

[file_input, pathname] = uigetfile({'*.sto', 'OpenSim Controls (excitation) Files (*.sto)'}, ...
                                   'Select the initial controls file','MultiSelect', 'off');
temp_i = importdata(strcat(pathname,file_input)); % import controls data (1st column is time)
old_time_i = temp_i.data(:,1);         % time 
old_data_i = temp_i.data(:,2:end);     % the controls
old_text_i = temp_i.textdata(7,2:end); % controls names


% Arrange the initial guess by nodes and controls

u0_temp = zeros(N,Ncontrols);

for j = 1:size(controls_all,1)
    for k = 1:size(old_data_i,2)
        if strcmp(old_text_i(k),controls_all(j)) == 1,
            % interpolate to the temporal grid
            u0_temp(:,j) = interp1(old_time_i,old_data_i(:,k),dc_time);
        end
    end
end

% Append the initial guess for the controls to the end of X0

for i = 1:Ncontrols,
    X0(Nstates*N + N*(i-1)+1 : Nstates*N + N*i,:) = u0_temp(:,i);
end

% Check: make sure both files (states and controls) are consistent

if (old_time_i(end) ~= old_time_s(end))
    disp('Time stamp of states and controls for initial guess do not match!');
    disp('Check the input files...');
    return
end

Pos_LB(1:Ncoord*N) = -0.2;       Pos_UB(1:Ncoord*N) = 0.2;

% Activation and fiber length are arranged in alternate manner

Vel_LB(1:Ncoord*N) = -2.0;       Vel_UB(1:Ncoord*N) = 2.0;

Act_LB(1:Nmuscles*N) = 0.011;    Act_UB(1:Nmuscles*N) = 0.999;
Fib_LB(1:Nmuscles*N) = 0.011;    Fib_UB(1:Nmuscles*N) = 0.999;

for i = 1:Nmuscles, % Activation and fiber length are arranged in alternate manner
    Act_Fib_LB(2*N*(i-1)+1:2*N*(i)) = [Act_LB(N*(i-1)+1:N*(i)) Fib_LB(N*(i-1)+1:N*(i))];
    Act_Fib_UB(2 *N*(i-1)+1:2*N*(i)) = [Act_UB(N*(i-1)+1:N*(i)) Fib_UB(N*(i-1)+1:N*(i))];
end

i = linspace(-0.5086,0.5086,150);
k = sin(i);

Con_LB(1:Nmuscles*N) = k;   

%r = sin(s);
%Con_UB1 = s;
%Con_UB2=s;
Con_UB=linspace(0.5,0.5,150);


lb = [Pos_LB Vel_LB Act_Fib_LB Con_LB]';
ub = [Pos_UB Vel_UB Act_Fib_UB Con_UB]';
 g = transpose(lb);
 b= transpose(ub);
 
 net = feedforwardnet(10,'traingd');
 %[net,tr] =train(net,g,b);%
 net = configure(net,g,b);
 y1=net(g);

 [net,tr] = train(net,g,b);
  plotperf(tr);
 s=net(g);
 a=reshape(b,1,[]);
  z= transpose(a);
 
 % Put Xopt back in 'row = time, column = states/controls' format

X_state_opt = zeros(N,Nstates); %pre-allocate size

for i = 1:Nstates,
    X_state_opt(:,i) = z(N*(i-1)+1:N*i,1);
end

X_controls_opt = zeros(N,Ncontrols); %pre-allocate size

for i = 1:Ncontrols,
    X_controls_opt(:,i) = z(Nstates*N + N*(i-1)+1:Nstates*N + N*i,1);
end

% Create data structure for the states error file

StatesData = struct();
StatesData.name = [char(osimModel.getName()), 'error_state', char(date)];
StatesData.nRows = size(dc_time, 1);
StatesData.nColumns = Nstates+1; %All the states + time
StatesData.inDegrees = false;
StatesData.labels = cell(1,StatesData.nColumns); 
StatesData.labels{1}= 'time';

for j = 2:1:StatesData.nColumns
   StatesData.labels{j} = char(states_all(j-1));
end
StatesData.data = [dc_time, X_state_opt];
writeOpenSimStatesFile(StatesData)

% Create data structure for the controls error file

ControlData = struct();
ControlData.name = [char(osimModel.getName()), 'error_control', char(date)];
ControlData.nRows = size(dc_time, 1);
ControlData.nColumns = Ncontrols+1; %All the controls + time
ControlData.inDegrees = false;
ControlData.labels = cell(1,ControlData.nColumns); 
ControlData.labels{1}= 'time';

for j = 2:1:ControlData.nColumns
   ControlData.labels{j} = char(controls_all(j-1));
end
ControlData.data = [dc_time, X_controls_opt];
writeOpenSimControlFile(ControlData)
