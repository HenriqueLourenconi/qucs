classdef synctrcircuit < qucstrans
    % synchonous circuit solver for use with the ode integrator routines
    % based on an interface to the qucs transient circuit solvers
    %
    %
    
    % Copyright Richard Crozier 2013
    
    properties (SetAccess = private, GetAccess = public)
        
        last_t;
        
    end
    
    methods
        
        function this = synctrcircuit(netlist)
            
            this = this@qucstrans(netlist);
            
        end
        
        function init(this, tstart)
            % initialises the circuit solver
            
            this.cppcall('init_sync', tstart);
            
            this.isinitialised = true;
            
            this.last_t = tstart;
            
        end
        
        function stepsolve(this, t)
            % solves a synchronous step at the supplied time
            %
            % Syntax
            %
            % stepsolve(this, t)
            %
            % Input
            %
            %   t - the time point at which the circuit simulation is to be
            %   solved
            %
            % Description
            %
            % This function causes the qucs solvers to attempt to solve the
            % circuit at the time t. The result is stored internally, and
            % can be retrieved with the getsolution method. The solution is
            % not added to the history of time steps used by the circuit
            % integrators. To add the solution to the history and move on
            % to the next step, the acceptstep_sync function must be used.
            %
            
            if ~isscalar(t)
                error('t must be a scalar.');
            end
            
            this.cppcall('stepsolve_sync', t);
            
        end
        
        function [sol, NM] = getsolution(this)
            % retrieves the last solution generated by the circuit solver
            %
            % Syntax
            %
            % [sol, NM] = getsolution(this)
            %
            % Description
            %
            % getsolution returns an vector containing the nodal voltages
            % and branch currents at the last solution generated by a call
            % to stepsolve acceptstep, or just initialisation. The number
            % of nodal voltages and branch currents respectively are
            % returned in the two element vector NM.
            %
            
            [sol, NM] = this.cppcall('getsolution');
            
        end
        
        function acceptstep(this, t)
            % recalculates and accepts a time step solution into the
            % internal history of solutions for the circuit integrators
            %
            % Syntax
            %
            % acceptstep(this, t)
            %
            % Input
            %
            %   t - the time point at which the circuit is to be
            %     recalculated and added to the solution history
            % 
            
            stepsolve(this, t);
            
            this.cppcall('acceptstep_sync');
            
            this.last_t = t;
            
        end
        
    end
    
    
    methods(Static)
        
        function status = odeoutputfcn(t, y, flag, qtr_sync)
            % static method for updating the synchronous circuit on each
            % time step of an ode solution
            %
            % Syntax
            %
            % status = odeoutputfcn(t, y, flag, qtr_sync)
            %
            % Description
            %
            % This function must be called as the output function, or
            % within the output function of an ode solution incorporating
            % the synchronous circuit solver. See the help for odeset.m for
            % more details on the ode solver OutputFcn option. 
            %
            % 
            
            status = 0;
            
            if isempty(flag)
                % accept the time step into the circuit solution history
                qtr_sync.acceptstep(t);
                
            elseif strcmp(flag, 'init')
                % initialise the circuit
                qtr_sync.init(t(1)+10*eps(t(1)));
                
            elseif strcmp(flag, 'done')
                
%                 qtr_sync.
                
            end
            
        end
        
        function dydt = odefcn(t, y, qtr_sync)
            % static method for determining the derivatives of the circuit
            % solution at each time step for use with the ode solvers
            %
            % Syntax
            %
            % status = odefcn(t, y, qtr_sync)
            %
            
            % test if the circuit is initialised, this is necessary,
            % because contrary to what the ode solver documentation states,
            % the solver function is called before the OutputFcn is ever
            % called with the 'init' flag in odearguments.m
            if qtr_sync.isinitialised
                % if the circuit is initialised start solving
                
                % attempt to solve the circuit at the current time step
                qtr_sync.stepsolve(t);
                
                % get the solution just calculated
                sol = qtr_sync.getsolution();
                
                % get the derivative w.r.t.
                dydt = (sol - y) / (t - qtr_sync.last_t);
                
            else
                dydt = zeros(length(y),1);
            end
            
        end
        
        
        
    end
    
end