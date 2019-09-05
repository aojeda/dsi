classdef RSBL < handle
    properties
        H
        Delta
        Blocks
        Ny
        Nx
        Ng
        defaultOptions = struct(...
            'maxIter',100,...       % Maximum number of iterations
            'verbose',true,...      % Produce per-iteration prints
            'maxTol',1e-1,...       % Maximum tolerance of logE change
            'bufferSize',100,...    % History buffer size
            'doPruning','true',...  % Enable the pruning stage
            'learnLambda','true',...% Enable the pruning stage
            'learnGammaF','true');  % Enable the pruning stage
    end
    properties(GetAccess=private)
        Hi
        HiHit
        CiHt
        Ci
        Ut
        s2
        Iy
        History
    end
    
    methods
        function obj = RSBL(H, Delta, blocks)
            obj.H = H;
            obj.Delta = Delta;
            obj.Blocks = blocks;
            [obj.Ny,obj.Nx] = size(H);
            obj.Ng = size(blocks,2);     
            obj.Hi    = cell(1,obj.Ng);
            obj.HiHit = zeros([obj.Ny, obj.Ny, obj.Ng]);
            obj.Ci    = sparse(obj.Nx^2,obj.Ng);
            obj.CiHt  = sparse(obj.Nx*obj.Ny,obj.Ng);
            ind1      = zeros(obj.Nx,obj.Nx);
            ind2      = zeros(obj.Nx,obj.Ny);
            obj.Iy    = speye(obj.Ny);
            obj.History = struct('lambda',nan(obj.defaultOptions.bufferSize,1),'gamma_F',nan(obj.defaultOptions.bufferSize,1),'logE',nan(obj.defaultOptions.bufferSize,1),'pointer',1);
            fprintf('Precomputing  ');
            for k=1:obj.Ng
                % Per-block square root precision matrix
                Di = Delta(blocks(:,k),blocks(:,k));
                
                % Per-block covariance matrix
                sqCi_k = inv(Di);
                ind1(blocks(:,k),blocks(:,k)) = 1;
                ind2(blocks(:,k),:) = 1;
                C_i = sqCi_k*sqCi_k';
                obj.Ci(ind1==1,k) = C_i(:);
                CiHt_k = C_i*obj.H(:,blocks(:,k))';
                obj.CiHt(ind2==1,k) = CiHt_k(:);
                ind2(blocks(:,k),:) = 0;
                ind1(blocks(:,k),blocks(:,k)) = 0;
                
                % Per-block standardized gain matrices
                obj.Hi{k} = obj.H(:,blocks(:,k))/Di;
                obj.HiHit(:,:,k) = obj.Hi{k}*obj.Hi{k}';
                
                if ~mod(k,round(obj.Ng/10)) 
                    fprintf('%i%% ',round(100*k/obj.Ng));
                end
            end
            fprintf('\n')
            
            % Unweighted prior covariance
            C = reshape(sum(obj.Ci,2),[obj.Nx, obj.Nx]);
                
            % Fix possible 0 diagonal elements
            dc = diag(C);
            dc(dc==0) = median(dc(dc~=0));
            C = C - diag(diag(C)) + diag(dc);
            
            % Compute svd
            sqC = chol(C);
            [U,s] = svd(obj.H*sqC,'econ');
            obj.s2 = diag(s).^2;
            obj.Ut = U';
            
            if ~exist('compileInvChol','file')
                addpath(fullfile(fileparts(which('RSBL')),'invChol'));
            end
            try
                invChol_mex(eye(4));
            catch
                compileInvChol;
            end
        end
        
        %%
        function [lambda, gamma, gamma_F, history] = learning(obj,Y, lambda0, gamma_F0, gamma0, options)
            if nargin < 6
                options = obj.defaultOptions;
            end
            if ~isfield(options,'doPruning')
                options.doPruning = true;
            end
            %[lambda0, gamma_F0] = initHyperparameters(obj, y);
            %UtY2 = (obj.Ut*Y).^2;
            %gamma_F0 = (UtY2'-lambda0)*(1./(obj.s2+eps));
            [lambda, gamma, gamma_F, history] = optimizeFullModel(obj,Y,lambda0, gamma_F0, options);
            if options.doPruning
                [gamma,history] = pruning(obj, Y, lambda, gamma, gamma0, history, options);
            end
        end
        
        %%
        function Sw = getSw(obj)
            Sw = reshape(obj.Ci*gamma,[obj.Nx obj.Nx]);
        end
        %%
        function [x, lambda, gamma_F, gamma, logE, history] = update(obj,y,lambda, gamma0, options)
            if nargin < 3
                lambda = [];
            end
            if nargin < 4
                gamma0 = eps;
            end
            if isempty(lambda)
                [lambda, gamma_F0] = initHyperparameters(obj, y);
            else
                [~, gamma_F0] = initHyperparameters(obj, y);
            end
            if nargin < 5
                options = obj.defaultOptions;
            end
            [lambda, gamma, gamma_F, history] = learning(obj, y, lambda, gamma_F0, gamma0, options);
            logE = calculateLogEvidence(obj,y,lambda,gamma);
            
            K = getK(obj, lambda, gamma);
            x = K*y;
        end
        %%
        function K = getK(obj, lambda, gamma)
            [~, iSy] = obj.calculateModelCov(lambda, gamma);
            SxHt = reshape(obj.CiHt*gamma,[obj.Nx obj.Ny]);
            K = SxHt*iSy;
        end
        %%
        function ypred = predict(obj, x)
            ypred = obj.H*x;
        end
        %%
        function logE = calculateLogEvidence(obj,y,lambda,gamma)
            [Sy, iSy] = calculateModelCov(obj, lambda, gamma);
            logE = (-1/2)*(y'*iSy*y + RSBL.logDet(Sy));
        end
    end
    methods(Access=private)
        %%
        function [lambda0, gamma0] = initHyperparameters(obj, Y)
            UtY2 = (obj.Ut*Y).^2;
            S = [obj.s2 obj.s2*0+1];
            phi = abs(mean((S'*S)\(S'*UtY2),2));
            gamma0  = phi(1);
            lambda0 = phi(2);
        end
        %%
        function [lambda, gamma, gamma_F, history] = optimizeFullModel(obj,Y,lambda0, gamma_F0, options)
            UtY2 = (obj.Ut*Y).^2;
            Nt = size(Y,2);
            gamma = ones(obj.Ng,1);
            lambda = lambda0;
            gamma_F = gamma_F0;
            gamma(:) = gamma_F;
            
            history = obj.History;
            history.lambda(1)  = lambda;
            history.gamma_F(1) = gamma_F;
            history.logE(1)    = calculateLogEvidence(obj,Y,lambda,gamma);
            
            for k=2:options.maxIter
                psi = gamma_F*obj.s2+lambda;
                psi2 = psi.^2;
                
                %lambda   = (lambda0+lambda *(sum(mean(bsxfun(@times,UtY2,     1./psi2),2)))/(eps+sum(     1./psi)))/2;
                gamma_F  = (gamma_F0+gamma_F*sum(mean(bsxfun(@times,UtY2,obj.s2./psi2),2))/(eps+sum(obj.s2./psi)))/2;
                %lambda   = lambda *sum(mean(bsxfun(@times,UtY2,     1./psi2),2))/(eps+sum(     1./psi));
                %gamma_F  = gamma_F*sum(mean(bsxfun(@times,UtY2,obj.s2./psi2),2))/(eps+sum(obj.s2./psi));
                gamma(:) = gamma_F;
                history.lambda(k) = lambda;
                history.gamma_F(k) = gamma_F;
                history.logE(k) = calculateLogEvidence(obj,Y,lambda,gamma);
                
                if options.verbose
                    fprintf('%i => diff(logE): %.4g   logE: %.5g   Lambda: %.4g   Gamma: %.4g\n',...
                        k,abs(diff(fliplr(history.logE(k-1:k)))),history.logE(k),lambda,gamma_F);
                end
                
                % Check convergence and exit condition
                if diff(history.logE(k-1:k)) < options.maxTol, break;end
            end
            history.pointer = k;
        end

        %%
        function [gamma, history] = pruning(obj,Y, lambda, gamma, gamma0, history, options)
            % Implements the \gamma-MAP SBL algorithm of David Wipf and B.
            % Rao.
            %
            %% -------------------------------------------------------------
            % For certain problems the following algorithm may be faster:
            %
            %--Precompute the following matrices in the constructor (once)--
            % obj.A = sparse([]);
            % for i=1:obj.Ng
            %     obj.A = blkdiag(obj.A,obj.Hi{i}');
            % end
            % obj.Ig = speye(obj.Ng);
            %
            %--Substitute this function with the following code-------------
            % Nt = size(Y,2); 
            % iSy_sq = sqrtm(iSy);
            % Yk = kron(obj.Ig,iSy_sq'*Y);
            % B = kron(obj.Ig,iSy_sq);
            % C = obj.A*B;
            % A = sum(reshape(sum((C*Yk).^2),Nt,obj.Ng))';
            % B = sum(reshape(sum(C.^2), obj.Ny, obj.Ng))';
            % gamma = (gamma0+gamma.*(A./B))/2;
            %% -----------------------------------------------------------------------------
            if numel(gamma0)==1
                gamma0 = gamma*0+gamma0;
            end
            A = 0*gamma;
            B = A;
            for k=1:options.maxIter
                [~, iSy] = obj.calculateModelCov(lambda,gamma);
                for i=1:obj.Ng
                    Hi_iSy = obj.Hi{i}'*iSy;
                    A(i) = norm(Hi_iSy*Y,'fro')^2;
                    B(i) = (abs(sum(sum((Hi_iSy)'.*obj.Hi{i}))));
                end
                gamma = 0.25*gamma0 + 0.75*gamma.*(A./B);
                history.pointer = history.pointer+1;
                history.logE(history.pointer) = calculateLogEvidence(obj,Y,lambda,gamma);
                if options.verbose
                    fprintf('%i => diff(logE): %.4g   logE: %.5g   Sum Gamma: %.4g\n',history.pointer,diff(...
                        history.logE(history.pointer-1:history.pointer)),history.logE(history.pointer),sum(nonzeros(gamma)));
                end
                if diff(history.logE(history.pointer-1:history.pointer)) < options.maxTol, break;end
            end
        end
        
        %%
        function [Sy, iSy] = calculateModelCov(obj,lambda,gamma)
            gHHt = sum(bsxfun(@times, obj.HiHit,permute(gamma,[3 2 1])),3);
            Sy = lambda*obj.Iy+gHHt;
            try
                iSy = invChol_mex(double(Sy));
            catch ME
                if ~strcmp(ME.identifier,'MATLAB:invChol_mex:iscomplex')    
                    warning(ME.message);
                end
                if strcmp(ME.identifier,'MATLAB:invChol_mex:dpotrf:notposdef')
                    warning('Possibly the data is rank deficient!')
                end
                [Utmp,S,Vtmp] = svd(Sy);
                stmp = real(diag(S));
                invS = 1./stmp;
                invS(isinf(invS)) = 0;
                iSy = Utmp*diag(invS)*Vtmp';
            end
        end
    end
    methods(Static)
        %%
        function log_d = logDet(S)
            log_d = log(det(S));     
            if isinf(log_d)
                e = eig(S);
                e(e<0) = eps;
                log_d = sum(log(e));
            end
        end
    end
end
