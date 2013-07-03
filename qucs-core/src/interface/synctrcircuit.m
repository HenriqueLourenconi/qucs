classdef synctrcircuit < qucstrans
    % synchonous circuit solver for use with the ode integrator routines
    % based on an interface to the qucs transient circuit solvers
    %
    %
    
    % Copyright Richard Crozier 2013
    
    properties
        
        
    end
    
    methods
        
        function this = synctrcircuit(netlist)
            
            this = this@qucstrans(netlist);
            
        end
        
        function init(this, tstart)
            this.cppcall('init_sync', tstart);
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
            
        end
        
        
    end
    
end