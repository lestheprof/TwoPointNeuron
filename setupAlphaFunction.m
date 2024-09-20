function alphafunction = setupAlphaFunction(timestep,tau)
            % how long does the alpha function need to be stored for?
            % according to Koch (p100), alpha function reaches 1 percent of
            % the peak at 7.64 * tpeak, and tpeak = tau. So we'll use 10 *
            % tau for a little safety.
            %
            % LSS last updated 16 July 2024
            %
            alphalength = round(10 * tau / timestep) ;
            alphafunction = zeros([1  alphalength]) ;
            for i=1:alphalength
                alphafunction(i) = (i*timestep/tau) * exp(1 - (i* timestep)/tau) ;
            end
            % normalise sum of alphafunction elements to 1
            alphafunction = alphafunction/sum(alphafunction) ; % normalise
        end