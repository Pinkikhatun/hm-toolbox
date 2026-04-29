classdef hsstangent
    %HSSTANGENT Summary of this class goes here
    %   Detailed explanation goes here

   properties
       %properties of the tangent vectors
        TB12
        TB21
        
        % Factorization of the upper triangular block as U12 * V12'
        %collinear component of U and V
        TU
        TV
        %%Orthogonal component with respect to U and V
        PU
        PV
        % Factorization of the lower triangular block as U21 * V21'
        TRl
        PRl
        TRr
        PRr
        TWl
        PWl
        TWr
        PWr
        
        % Dense version of the matrix, if the size is smaller than the
        % minimum allowed block size.
        TD
        
        
        % top and bottom blocks of the matrix.
        B12
        B21
        

        %Part related HSS matrix
        % Factorization of the upper triangular block as U12 * V12'
        U
        V
        
        % Factorization of the lower triangular block as U21 * V21'
        Rl
        Rr
        Wl
        Wr
        
        % Size of the matrix
        ml
        nl
        mr
        nr
        
        % Dense version of the matrix, if the size is smaller than the
        % minimum allowed block size.
        D
        
        topnode
        leafnode
        
        TA11
        TA22
        
    end

    methods
        function obj = hsstangent(varargin)
            if nargin == 0
                return;
            end    
            if ischar(varargin{1})
            switch varargin{1}
                case 'randn'
                H=varargin{2};
                r=varargin{3};
                obj= hsstangent();
                obj= create_randvec(H, obj,r);
                otherwise
                    error('unsupported constructor mode');
            end
            else
                % case where the vector is specified by a cell array
                obj= hsstangent();
                obj=create_randvec_from_cell(varargin{1}, varargin{2}, obj);
            end

        end
    end
end