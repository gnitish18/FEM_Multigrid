%% Function to display the animation of Beam Deflection
%   Input Parameters:
%   * U_fem - FEM solution deflection matrix
%   * U_ana - Analytical solution deflrction matrix
%   * X     - Nodal locations along beam length
%%

function Animate(U_fem, U_ana, X)
 
    % Define number of frames
    framemax = 90;
    M = moviein(framemax);
    
    % Set the position of the animation window
    set(gcf,'Position',[100 100 640 480]);

    for i = 1:framemax

        % Interpolate the beam deflections framewise
        u_fem = U_fem*i/framemax;
        u_ana = U_ana*i/framemax;
        
        % Plot the deflection alog the length of the beam
        figure(4)
        plot(X,u_fem,'k','LineWidth',5);
        hold on
        scatter(X,u_ana,'y','LineWidth',1.5);
        axis([0 5 -4000 1500])
        title('Deflection of Cantilever Beam','fontsize',18)
        ylabel('Deflection of Beam (units)') 
        xlabel('Length of Beam (units)') 
        legend('FEM','Analytical')
        legend('boxoff')
        hold off
        
        % Add frames to list
        M(:,i) = getframe(gcf);
        if i == framemax
            pause(1.5)
        end

    end
    
    % Reset graphics properties
    clf reset
    axis off
    
    % Close animation
    close(figure(4))
    
end