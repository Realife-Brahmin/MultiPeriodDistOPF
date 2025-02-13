% clc;
% clear all;
function [V, x, Table,Volttable]= singlephaselin()
load linedata.txt        % load the line data, it has R and X value for the edges
branch=linedata;
load powerdata10.txt;    % load the power data, node number, P load, Q load, Capacitor, DG
powerdata = powerdata10;
% load powerdata50.txt;    % load the power data, node number, P load, Q load, Capacitor, DG
% powerdata = powerdata50;
mult = 1;                                % load mupltiplier
mult1 = 1;                               % PV mupltiplier

fb = branch(:,1);                        % from bus
tb = branch(:,2);                        % to bus
G = graph(fb,tb);                        % generate the tree of the network,
nb = size((powerdata),1);               % total number of nodes 
bKVA = 100000;                             % base kVA

resitance = ((branch(:,3)));          % resistance of the edge
reactance = ((branch(:,4)));          % reactance of the edge

PL = powerdata(:,2)/bKVA;              % rated active power of node
QL = powerdata(:,3)/bKVA;              % rated reactive power of node
QC = zeros(nb,1);
PG = powerdata(:,4)/bKVA;              % rated active power of node


PL = PL.*mult;                          % rated active power multipiled with load multiplier                     
QL = QL.*mult;                           % rated active power multipiled with load multiplier                     
PG = PG.*mult1;                          % rated active power multipiled with load multiplier  

T=dfsearch(G,1,'edgetonew');                % finding Parent child combination using dfsearch

%% calculating R and X matrix
R = zeros(nb);
X  = zeros(nb);
for i = 1:(nb-1)
    R(fb(i), tb(i)) = resitance(i);
    R(tb(i) ,fb(i)) = R(fb(i), tb(i)) ;
    X(fb(i), tb(i))= reactance(i);
    X(tb(i) ,fb(i)) = X(fb(i), tb(i)) ;    
end

%% Finding number of unknowns

 Ap = 1:nb-1;                          % defining the variabes for P
 Aq = nb:2*(nb-1);                      % defining the variabes for Q
 Av = 2*(nb):3*(nb-1)+1;                % defining the variabes for V

Table = [T(:,1) T(:,2) Ap'  Aq'   Av'];     % creating Table for variabes P, Q , V

Da = 2*(nb)-1:3*(nb-1)+1;               
Volttable = Da';                        % voltage variables including substation


for i = 1:(nb)                       % number of child of a node 
tc(i)=size(find((i)==T(:,1)),1) ;
row = find(i == T(:,1));
end
tc;


%% Formation of equality constarints
% Aeq = zeros(90,90);
% beq = zeros(90,1);
Vs= 1.0;                                % substation voltage
 for i =2:nb
    k = tc(i) ;                       %% total # of child
    row = find(i == T(:,1));
 
  if isempty(row)                       % if there is no child for a node
    Parent = find(i == T(:,2));         % find parent of a node
    
    if isempty(Parent)  
        Aeq = Aeq;

    else    
        Poc = find(T(Parent,1) == T(:,1));
        Aeq(Parent,Table(Parent,3))= 1;                                               % P12 = PL2
        Aeq(Parent+(nb-1),Table(Parent,4))= 1;                                        % Q12 = QL2-QC2
        Aeq(Parent+2*(nb-1),Table(Parent,5))= 1;                                      % V2 -V1 + 2rP+2xP =0                                           
        Aeq(Parent+2*(nb-1),Volttable(Poc(1)))= -1;
        Aeq(Parent+2*(nb-1),Table(Parent,3))= 2*(R(T((Parent),1),T((Parent),2)));
        Aeq(Parent+2*(nb-1),Table(Parent,4))= 2*(X(T((Parent),1),T((Parent),2)));
        beq(Parent)= 1*PL(Table(Parent,2))-PG(Table(Parent,2));
        beq(Parent+(nb-1)) = 1*QL(Table(Parent,2))-QC(Table(Parent,2));
   
    end
 
  else                                                      % if the node has child nodes 
    Parent = find(i == T(:,2));
    
    if isempty(Parent)  
        Aeq = Aeq;
    
    else    
        Poc = find(T(Parent,1) == T(:,1));
        Aeq(Parent,Table(Parent,3))= 1;                             % P12- summation(P23) = PL2 
        Aeq(Parent,Table(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(Parent, Table(row(j+1),3)) =   - 1;
            end
            
        Aeq(Parent+(nb-1),Table(Parent,4))= 1;                      % Q12- summation(Q23) = QL2-QC2 
        Aeq(Parent+(nb-1),(nb-1)+ Table(row(1),3))= -1;
            for j = 1:length(row)-1
                Aeq(Parent+(nb-1),(nb-1)+Table(row(j+1),3)) =  -  1;
            end
   
       Aeq(Parent+2*(nb-1),Table(Parent,5))= 1;                                 % V2 -V1 + 2rP+2xP =0  
       Aeq(Parent+2*(nb-1),Volttable(Poc(1)))= -1;
       Aeq(Parent+2*(nb-1),Table(Parent,3))= 2*(R(T((Parent),1),T((Parent),2)));
       Aeq(Parent+2*(nb-1),Table(Parent,4))= 2*(X(T((Parent),1),T((Parent),2)));
       beq(Parent)= 1*PL(Table(Parent,2))-PG(Table(Parent,2));
       beq(Parent+(nb-1)) =  1*QL(Table(Parent,2))-QC(Table(Parent,2));
     end
  end
  
 Aeq(3*(nb-1)+1,Volttable(1)) = 1;                                  % equation for substaion voltage
 beq(3*(nb-1)+1) = Vs;
 
 end
Aeq;

nvar = size(Aeq,2);                                             % totale number of variables
 

%% formation of objective function

f = zeros(nvar,1);

%% limits defination
lb(1:1458,1)= -(1500*ones(1458,1));
lb(1459:2188,1)= ((0.9^2)*ones(730,1));

ub(1:1458,1)= (1500*ones(1458,1));
ub(1459:2188,1)= ((1.1)*ones(730,1));
%%
optimoptions(@intlinprog,'Display','iter');
[x,fval,exitflag,output] = intlinprog(f,[],[],[],Aeq,beq,[],[]);

 V = sqrt(x(1459:2188));                  % voltage 
 Ptotal = x(1)*bKVA;                         % substaion active power
 Qtotal = x(730)*bKVA;                       % substaion reactive power

Vlin = 1;        %% substaion voltage
for i = 2:730
indx = find(i==Table(:,2));
Vnon1 =  sqrt(x(Table(indx,5)));
Vlin = [Vlin Vnon1];
end