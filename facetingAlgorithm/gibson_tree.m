function facets = gibson_tree(i, j, fracs, comps, facets, flb, bin, fu, ac, outcome)


%%

%% Begin determining bin
if fracs(1) > 0.5 % set facet to bin with 50% or more area
    facets(i,j) = comps(1);
    if outcome; disp('test 2'); end
    
else
    if sum(fracs(1:2))>=0.5 % if bin 1&2 contain 50% of area
        flb = NaN(3,1); % flatness
        for m = 1:3
            % flb(m) = percent of area in each direction that is flat
            flb(m) = sum(fu(bin==comps(m)))/numel(fu(bin==comps(m))); %%%%%%THIS GOT EDITED numel(bin)
        end
        % true if first direct is adjacent to second or third
        adj1_23 = orientation_is_adjacent(comps,1,[2 3]);


        %% check if (ANY DIRECT IS 50% FLAT) or (DIFF B/W 1&2 IS >
        %% 0.2) or (adjacent facets)
        if any(flb)>=0.5 || abs(ac(1)-ac(2))>=0.2 || adj1_23 % check 2------------
            facets(i,j) = comps(1);
            if outcome; disp('rule 1'); end


        else % check 2 no, goto check 3
            % check if second facet is adjacent to third facet
            adj2_3 = orientation_is_adjacent(comps,2,3);

            % if not 1 adj to 2or3 and not 2 adj to 3
            if ~adj1_23 && ~adj2_3 % check 3 -------------------
                facets(i,j) = comps(1);
                if outcome; disp('rule 2'); end


            else % check 3 no, goto check 4
                adj1_4 = orientation_is_adjacent(comps,1,4);
                adj2_3 = orientation_is_adjacent(comps,2,3);

                %% if 1 is adj to 4 && 2 adj to 3
                if adj1_4 && adj2_3 % check 4 -------------------


                    % if diff b/w 4 & 3 is < 0.1
                    if abs(fracs(4)-fracs(3))<=0.1 % check 5 -------------------
                        % if 1&4 is larger than 2&3, set as 1
                        if fracs(1)+fracs(4) > fracs(2)+fracs(3) % check 6 -------------
                            facets(i,j) = comps(1);
                            if outcome; disp('rule 3'); end
                        else % if 2&3 is larger than 1&4, set as 2
                            facets(i,j) = comps(2);
                            if outcome;  disp('rule 4'); end
                        end

                    else % if diff b/t 4 & 3 is > 0.1
                        % if sum of 2&3 is larger than 1, set as 2
                        if sum(fracs(2:3)) > fracs(1) % check 8 ----------------
                            facets(i,j) = comps(2);
                            if outcome; disp('rule 5'); end
                        else % if sum of 2&3 is smaller than 1, set as 1
                            facets(i,j) = comps(1);
                            if outcome; disp('rule 6'); end
                        end
                    end


                %% if not (1 is adj to 4 && 2 is adj to 3)
                else % check 4 no, goto check 7
                    % if 2 is adj 3
                    if orientation_is_adjacent(comps,2,3) % check 7 --------------
                        % if 2 + 3 is greater than 1, set as 2
                        if fracs(2)+fracs(3)>fracs(1) % check 9 ----------------
                            facets(i,j) = comps(2);
                            if outcome; disp('rule 7'); end
                        else % if 2 + 3 is less than 1, set as 1
                            facets(i,j) = comps(1);
                            if outcome; disp('rule 8'); end
                        end

                    else % check 7 no --> error
                        error('logic error, not possible');
                    end
                end

            end
        end



    %% if bin 1 & 2 doesn't contain 50% of the area
    else 
        if isempty(flb)  % create flb
            flb = NaN(3,1);
            for m = 1:3
                flb(m) = sum(fu(bin==comps(m)))/numel(fu(bin==comps(m))); %%%%%%%% THIS GOT EDITED numel(bin);
            end
        end

        %% if 1&2 are adj && 1 is less than 50% flat && 2 is less
        %% than 50% flat, set as 1
        if orientation_is_adjacent(comps,1,2) && flb(1)<0.5 && flb(2)<0.5 % check 10 --------------
            % check 10 yes --> rule 10
            facets(i,j) = comps(1);
            if outcome; disp('rule 10'); end

        %%    
        else
            %% if 1 is less than 50% flat && 2 is less than 50% flat
            %% && 1&2 are NOT adj
            if flb(1)<0.5 && flb(2)<0.5 && ~orientation_is_adjacent(comps,1,2) % check 11 ------------------
                % check 11 yes, goto check 12

                if orientation_sum(comps,fracs,1) > orientation_sum(comps,fracs,2) % check 12 ------------
                    facets(i,j) = comps(1);
                    if outcome; disp('rule 11'); end
                else 
                    facets(i,j) = comps(2);
                    if outcome; disp('rule 12'); end
                end

            %%   
            else % check 11 no, goto check 13
                if flb(1)>=0.5 || flb(2)>=0.5 % check 13
                    % check 13 yes, goto check 14
                    if flb(1)>0.5 % check 14

                        if orientation_sum(comps,fracs,2) > fracs(1) % check 15 ---------
                            facets(i,j) = comps(2);
                            if outcome; disp('rule 13'); end
                        else % check 15 no --> rule 14
                            facets(i,j) = comps(1);
                            if outcome; disp('rule 14'); end
                        end

                    else % check 14 no, goto check 16
                        if fracs(1)>=fracs(2) % check 16
                            facets(i,j) = comps(1);
                            if outcome; disp('rule 15'); end
                        else % check 16 no --> rule 16
                            facets(i,j) = comps(2);
                            if outcome; disp('rule 16'); end
                        end

                    end
                else % check 13 no --> error
                    error('logic error, not possible');
                end
            end
        end                
    end  
end   
            