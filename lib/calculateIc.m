function Ic = calculateIc(Dp, Dq, Do, Nt, Area, expDo, Co, Lt, DD, expDq, m, expDp)
    % CALCULATEIC Compute the Physical Optics area integral Ic for a triangle.
    %
    %   Handles 4 special cases based on phase term magnitudes:
    %     Case 1: |Dp| < Lt, |Dq| >= Lt
    %     Case 2: |Dp| < Lt, |Dq| < Lt
    %     Case 3: |Dp| >= Lt, |Dq| < Lt
    %     Case 4: |Dp| >= Lt, |Dq| >= Lt, |DD| < Lt
    %     General case otherwise
    %
    %   Input:
    %       Dp, Dq, Do - phase terms
    %       Nt         - number of Taylor series terms
    %       Area       - triangle area array
    %       expDo      - exp(1i * Do)
    %       Co         - wave amplitude (1.0)
    %       Lt         - Taylor series threshold
    %       DD         - Dq - Dp
    %       expDq      - exp(1i * Dq)
    %       m          - triangle index
    %       expDp      - exp(1i * Dp)
    %   Output:
    %       Ic - complex area integral result

    if abs(Dp) < Lt && abs(Dq) >= Lt
        % Special case 1: |Dp| small, |Dq| large
        sic = 0.0;
        for n = 0:Nt
            sic = sic + (1i * Dp)^n / factorial(n) * (-Co / (n + 1) + expDq * (Co * GFunc(n, -Dq)));
        end
        Ic = sic * 2 * Area(m) * expDo / (1i * Dq);

    elseif abs(Dp) < Lt && abs(Dq) < Lt
        % Special case 2: both |Dp| and |Dq| small
        sic = 0.0;
        for n = 0:Nt
            for nn = 0:(Nt - 1)
                sic = sic + (1i * Dp)^n * (1i * Dq)^nn / factorial(nn + n + 2) * Co;
            end
        end
        Ic = sic * 2 * Area(m) * expDo;

    elseif abs(Dp) >= Lt && abs(Dq) < Lt
        % Special case 3: |Dp| large, |Dq| small
        sic = 0.0;
        for n = 0:Nt
            sic = sic + (1i * Dq)^n / factorial(n) * Co * GFunc(n + 1, -Dp) / (n + 1);
        end
        Ic = sic * 2 * Area(m) * expDo * expDp;

    elseif abs(Dp) >= Lt && abs(Dq) >= Lt && abs(DD) < Lt
        % Special case 4: both large but |DD| small
        sic = 0.0;
        for n = 0:Nt
            sic = sic + (1i * DD)^n / factorial(n) * (-Co * GFunc(n, Dq) + expDq * Co / (n + 1));
        end
        Ic = sic * 2 * Area(m) * expDo / (1i * Dq);

    else
        % General case
        Ic = 2 * Area(m) * expDo * (expDp * Co / (Dp * DD) - expDq * Co / (Dq * DD) - Co / (Dp * Dq));
    end
end
