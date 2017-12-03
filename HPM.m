% @author: Maziar Raissi

classdef HPM
    properties
        dt % time step size
        X1, U1 % data at time n
        X0, U0 % data at time n-1
        hyp % hyper-parameters
        NLML % negative log marginal likelihood
    end
    
    methods
        function obj = HPM(X1, U1, X0, U0, dt, hyp)
            obj.X1 = X1;
            obj.U1 = U1;
            
            obj.X0 = X0;
            obj.U0 = U0;
            
            obj.dt = dt;
            
            obj.hyp = hyp;
            
            fprintf('Total number of parameters: %d\n', length(obj.hyp));
        end
        
        function [NLML, D_NLML] = likelihood(obj, hyp)
                        
            y = [obj.U1(:); obj.U0(:)];
            
            X1_ = obj.X1;
            X0_ = obj.X0; U0_ = obj.U0;
            dt_ = obj.dt;
            
            sigma = exp(hyp(end));
            hyp_ = hyp(1:end-1);
            
            N = size(y,1);
            
            K11 = k11(X1_, X1_, hyp_, 0);
            K10 = k10(X1_, X0_, hyp_, U0_, dt_, 0);
            K00 = k00(X0_, X0_, hyp_, U0_, U0_, dt_, 0);
            
            K = [K11  K10;
                 K10' K00];
            
            % Cholesky factorisation
            L = jit_chol(K + sigma*eye(N));
            
            alpha = L'\(L\y);
            NLML = 0.5*y'*alpha + sum(log(diag(L))) + log(2*pi)*N/2;
            
            D_NLML = 0*hyp;
            Q =  L'\(L\eye(N)) - alpha*alpha';
            for i=1:length(hyp_)
                DK11 = k11(X1_, X1_, hyp_, i); 
                DK10 = k10(X1_, X0_, hyp_, U0_, dt_, i);
                DK00 = k00(X0_, X0_, hyp_, U0_, U0_, dt_, i);
                
                DK = [DK11  DK10;
                      DK10' DK00];
                
                D_NLML(i) = sum(sum(Q.*DK))/2;
            end
            
            D_NLML(end) = sigma*trace(Q)/2;
        end
        
        function obj = train(obj, n_iter)
                        
            [obj.hyp,~,~] = minimize(obj.hyp, @obj.likelihood, -n_iter);
            obj.NLML = obj.likelihood(obj.hyp);            
            
        end
        
        function [u1_star_mean, u1_star_var] = predict(obj, X1_star)
                        
            y = [obj.U1(:); obj.U0(:)];
            
            X1_ = obj.X1;
            X0_ = obj.X0; U0_ = obj.U0;
            dt_ = obj.dt;
            
            sigma = exp(obj.hyp(end));
            hyp_ = obj.hyp(1:end-1);
            
            N = size(y,1);
            
            K11 = k11(X1_, X1_, hyp_, 0);
            K10 = k10(X1_, X0_, hyp_, U0_, dt_, 0);
            K00 = k00(X0_, X0_, hyp_, U0_, U0_, dt_, 0);
            
            K = [K11  K10;
                 K10' K00];
            
            % Cholesky factorisation
            L = jit_chol(K + sigma*eye(N));
            
            K11 = k11(X1_star, X1_, hyp_, 0);
            K10 = k10(X1_star, X0_, hyp_, U0_, dt_, 0);
            
            psi = [K11 K10];
            
            u1_star_mean = psi*(L'\(L\y));
            
            alpha = (L'\(L\psi'));
            
            u1_star_var = k11(X1_star, X1_star, hyp_, 0) - psi*alpha;
        end
    end
end
