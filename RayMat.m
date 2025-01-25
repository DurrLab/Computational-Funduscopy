classdef RayMat
    % RayMat: A custom class for performing ray matrix calculations
    % for simple optical systems.
    %
    % This class encapsulates methods to represent, construct, and
    % manipulate optical systems using 2x2 ray matrices (Also called ABCD
    % matrices). Optical systems are defined as a linear sequence of
    % ordered optical elements.
    %
    % Properties:
    %   elements  - A cell array of 2x2 matrices representing individual
    %              optical elements in the system (ordered from input to
    %              output).
    %   rayMatrix - A single 2x2 matrix representing the overall transfer
    %               function of the optical system.
    %
    % Methods:
    %   RayMat(elements)      - Constructor to initialize the optical system.
    %   addElements(elements) - Add new optical elements to the system.
    %   updateRayMatrix()     - Update the overall system matrix (rayMatrix) based
    %                           on the current list of elements.
    %
    % Static Methods:
    %   thinlens(f)           - Returns the matrix for a thin lens of focal
    %                           length f.
    %   freespace(d)          - Returns the matrix for free-space propagation
    %                           over distance d.
    %   flatRefraction(n1,n2) - Returns the matrix for refraction at a flat
    %                           interface between media with refractive indices
    %                           n1 and n2.
    %   curvedRefraction(n1,n2,r) - Returns the matrix for refraction at a
    %                           curved interface with radius of curvature r.
    %   flatMirror()          - Returns the matrix for reflection off a flat
    %                           mirror.
    %   thickLens(n1,n2,r1,r2,t) - Returns the matrix for a thick lens defined
    %                           by its radii of curvature (r1, r2), refractive
    %                           indices (n1, n2), and thickness (t).
    %   computeRayMatrix(elements) - Combines a sequence of 2x2 element matrices
    %                                into a single system matrix.
    %
    % Author: Corey Simmerer
    % Date: 1/8/2025
    %

    properties
        elements; % A cell array consisting of optical elements in first to last order.
        rayMatrix;      % A single 2x2 matrix representing the entire system.
    end

    methods
        % Constructor
        function obj = RayMat(elements)
            %RayMat Construct an instance of this class
            %   elements: A cell array of 2x2 matrices representing optical
            %   systems
            obj.elements = elements;
            % Update the system ray matrix
            obj = updateRayMatrix(obj);
        end

        function obj = addElements(obj,elements)
            %addElement appends an optical element to the end of the system
            %(2x2)
            %   elements: A cell array of 2x2 matrices representing optical
            %   systems
            L = length(elements);
            obj.elements{end+1:end+L} = elements;
            % Update the system ray matrix
            obj = updateRayMatrix(obj);
        end
        
        function obj = updateRayMatrix(obj)
            obj.rayMatrix = RayMat.computeRayMatrix(obj.elements);
        end
        
    end
    

    % Static methods
    methods(Access=public, Static=true)
        %% Basic ray matrices
        function M = thinlens(f)
            M = [1 0; -1/f 1];
        end

        function M = freespace(d)
            M = [1 d; 0 1];
        end

        function M = flatRefraction(n1,n2)
            M = [1, 0; 0 n1/n2];
        end

        function M = curvedRefraction(n1, n2, r)
            M = [1 0; (n1-n2)/(r*n2), n1/n2];
        end

        function M = flatMirror()
            M = eye(2);
        end

        function M = thickLens(n1, n2, r1, r2, t)
            M1 = curvedRefraction(n1,n2,r1);
            M2 = freespace(t);
            M3 = curvedRefraction(n2,n1,r2);
            M = M3*M2*M1;
        end
        
        %% Useful functions 
        function M = computeRayMatrix(elements)
            % computeRayMatrix Convert list of optical element matrices
            % into single matrix encapsulating the system.
            %   elements: A cell array of 2x2 matrices representing optical
            %   systems in forward order

            M = eye(size(elements{end}));
            for k = length(elements):-1:1
                M = M * elements{k};
            end
        end
    end
end