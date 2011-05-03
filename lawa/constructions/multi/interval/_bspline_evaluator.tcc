namespace lawa {

//--- linear evaluators --------------------------------------------------------

template <typename T>
T
_linear_bspline_inner_evaluator0(T x, unsigned short deriv)
{
    double value = 0.0;
    if(deriv == 0){
        if(0 <= x && x < 0.25){
            value = 2.171240593367237661667 - 11.00419101450331474149* x;
        } else if (-0.75 <= x &&  x < -0.5){
            value = -2.458152258280902698695 - 3.441413161593263778173* x;
        } else if (0.5 <= x && x < 0.75){
            value = 0.484297173401595733543 - 0.678016042762234026960* x;
        } else if (0.75 <= x && x < 1.0){
            value = -0.096859434680319146709 + 0.096859434680319146709* x;
        } else if (-1.0 <= x && x < -0.75){
            value = 0.491630451656180539739 + 0.491630451656180539739* x;
        } else if ( 0.25 <= x && x < 0.5){
            value = -1.30490347253766076747 + 2.90038524911627897507 * x;
        } else if (-0.25 <= x && x < 0){
            value = 2.171240593367237661667 + 5.70778203747481756346* x;
        } else if (-0.5 <= x && x < -0.25){
            value = 2.226035845481337351211 + 5.92696304593121632164 * x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1){
        if(0 <= x && x < 0.25){
            value = - 11.00419101450331474149 ;
        } else if (-0.75 <= x &&  x < -0.5){
            value = - 3.441413161593263778173 ;
        } else if (0.5 <= x && x < 0.75){
            value = - 0.678016042762234026960 ;
        } else if (0.75 <= x && x < 1.0){
            value = 0.096859434680319146709 ;
        } else if (-1.0 <= x && x < -0.75){
            value = 0.491630451656180539739 ;
        } else if ( 0.25 <= x && x < 0.5){
            value = 2.90038524911627897507 ;
        } else if (-0.25 <= x && x < 0){
            value = 5.70778203747481756346 ;
        } else if (-0.5 <= x && x < -0.25){
            value = 5.92696304593121632164 ;
        } else{
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_bspline_inner_evaluator1(T x, unsigned short deriv)
{
    double value = 0.0;
    if(deriv == 0){
        if(0.5 <= x && x < 1.0){
            value = 3.464101615137754587055 - 3.464101615137754587055*x;
        } else if(0.0 <= x && x < 0.5){
            value = 3.464101615137754587055 * x ;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1){
        if(0.5 <= x && x < 1.0){
            value =  - 3.464101615137754587055;
        } else if(0.0 <= x && x < 0.5){
            value = 3.464101615137754587055;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_bspline_inner_evaluator2(T x, unsigned short deriv)
{
    double value = 0.0;
    if(deriv == 0){
        if(0.25 <= x && x < 0.5){
            value = 5.89923827193725729450 - 14.13397337635908350175 *x;
        } else if (0.0<= x && x < 0.25){
            value = 9.46297971138994567625 *x;
        } else if (0.75 <= x && x < 1.0) {
            value = -1.677990269774715967092 + 1.677990269774715967092 * x;
        } else if (0.5 <= x && x < 0.75) {
            value = -2.664250113839495385577 + 2.99300339519442185840 * x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1){
        if(0.25 <= x && x < 0.5){
            value = - 14.13397337635908350175;
        } else if (0.0<= x && x < 0.25){
            value = 9.46297971138994567625 ;
        } else if (0.75 <= x && x < 1.0) {
            value =  1.677990269774715967092;
        } else if (0.5 <= x && x < 0.75) {
            value =  2.99300339519442185840 ;
        } else {
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_bspline_left_evaluator0(T x, unsigned short deriv)
{
    double value = 0.0;
    if(deriv == 0){
        if(0 <= x && x < 0.25){
            value = 2.171240593367237661667 - 11.00419101450331474149* x;
        } else if (0.5 <= x && x < 0.75){
            value = 0.484297173401595733543 - 0.678016042762234026960* x;
        } else if (0.75 <= x && x < 1.0){
            value = -0.096859434680319146709 + 0.096859434680319146709* x;
        } else if ( 0.25 <= x && x < 0.5){
            value = -1.30490347253766076747 + 2.90038524911627897507 * x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1){
        if(0 <= x && x < 0.25){
            value = - 11.00419101450331474149 ;
        } else if (0.5 <= x && x < 0.75){
            value = - 0.678016042762234026960 ;
        } else if (0.75 <= x && x < 1.0){
            value = 0.096859434680319146709 ;
        } else if ( 0.25 <= x && x < 0.5){
            value = 2.90038524911627897507 ;
        } else{
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

template <typename T>
T
_linear_bspline_right_evaluator0(T x, unsigned short deriv)
{
    double value = 0.0;
    if(deriv == 0){
        if (-0.75 <= x &&  x < -0.5){
            value = -2.458152258280902698695 - 3.441413161593263778173* x;
        } else if (-1.0 <= x && x < -0.75){
            value = 0.491630451656180539739 + 0.491630451656180539739* x;
        } else if (-0.25 <= x && x < 0){
            value = 2.171240593367237661667 + 5.70778203747481756346* x;
        } else if (-0.5 <= x && x < -0.25){
            value = 2.226035845481337351211 + 5.92696304593121632164 * x;
        } else {
            value = 0.0;
        }
    } else if (deriv == 1){
        if (-0.75 <= x &&  x < -0.5){
            value = - 3.441413161593263778173 ;
        } else if (-1.0 <= x && x < -0.75){
            value = 0.491630451656180539739 ;
        } else if (-0.25 <= x && x < 0){
            value = 5.70778203747481756346 ;
        } else if (-0.5 <= x && x < -0.25){
            value = 5.92696304593121632164 ;
        } else{
            value = 0.0;
        }
    } else {
        value = 0.0;
    }
    return value;
}

} // namespace lawa
