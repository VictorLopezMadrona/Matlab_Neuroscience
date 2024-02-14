function newFigure_subplot(h1,~)
%% Function to act on subplot click
%         Mouse click: Plots the selected subplot to a new figure
%  Ctrl + Mouse click: Delete subplot
    switch get(gcf,'SelectionType')
        case 'normal'
            F = figure();
            copyobj(h1.Children,gca(F));
            % Copy the selected subplot title
            tmp = get(h1,'title'); tmp = tmp.String;
            % Set title to the new figure
            title(gca(F), tmp);
            XLim = get(h1,'XLim');
            YLim = get(h1,'YLim');
            set(gca, 'XLim',XLim)
            set(gca, 'YLim',YLim)
        case 'alt'
            delete(h1);
            
            if false
                % Get the list of remaining channels (stored as titles)
                fig=flipud(findobj(gcf,'type','axes'));
                % Store titles to a cell array
                subplot_titles = cell(length(fig),1);
                for i=1:length(fig)
                    tmp = get(fig(i),'title');
                    subplot_titles(i) = tmp.String;
                end
                % Coma delimited channels
                %sprintf('%s,' , subplot_titles{:})
            end
    end
end