from flask import Flask, request, render_template
from sympy import *

app = Flask(__name__)


# Các hàm tính toán cơ bản
def f_fomat(f):
    f_ln = f.replace(log, lambda arg: Function('ln')(arg))
    f_latex = latex(f_ln)
    return f_latex


def f_fomat1(f, subs_dict=None):
    f_str = str(sympify(f)).replace('log', 'ln')
    f_sym = sympify(f_str)
    if subs_dict:
        f_sym = f_sym.subs(subs_dict)
    return latex(f_sym)

def froundf(value):
    return f"{float(value):.6f}".rstrip('0').rstrip('.')

def calc(f, a):
    f_calc = sympify(f).subs(symbols('x'), a)
    return f_calc

def tinhtong(f, a, n):
    s = 0
    for i in range(0, n, 1):
        s += calc(f, a[i])
    return round(s, 6)

def calc_two_var(f, a, b):
    x, y = symbols('x y')
    f_calc = sympify(f).subs({x: a, y: b})
    return f_calc

def solve_h2(f1, f2):
    a, b = symbols('a b')
    result = solve([sympify(f1), sympify(f2)], (a, b))
    return result

def in_mang(s, a):
    result = f"{s} "
    for element in a:
        result += f"{element:<9}"
    return result

def check(arr, sll):
    for i in range(sll):
        for j in range(i + 1, sll):
            if float(arr[i]) == float(arr[j]):
                return True
    return False

@app.route("/", methods=["GET", "POST"])
def home():
    result_simpson = None
    result_newton = None
    result_regression = None 
    result_euler = None
    table_rows = []
    regression_table = [] 
    active_tab = "Simpson"
    error = None

    if request.method == "POST":
        if 'simpson' in request.form:
            # [Logic Simpson giữ nguyên]
            f_input = request.form.get("f_input")
            a = float(request.form.get("a"))
            b = float(request.form.get("b"))
            n = int(request.form.get("n"))

            f = sympify(f_input)
            h = (b - a) / float(n)
            formatted_f = f_fomat(f)
            h_latex = r"\frac{b - a}{n}"
            integral_expr = r"I = \frac{h}{3}\left[ \text{(đầu + cuối)} + 2\Sigma_{\text{chẵn}} + 4\Sigma_{\text{lẻ}} \right]"



            s_chan = 0
            s_le = 0
            j = 0
            i = a

            while j < (n + 1):
                row = {"j": j, "x_i": round(i, 6), "even": "", "odd": ""}
                if j % 2 == 0:
                    row["even"] = round(float(calc(f, i)), 6)
                    if i > a and j < n:
                        s_chan += calc(f, i)
                if j % 2 != 0:
                    row["odd"] = round(float(calc(f, i)), 6)
                    if i < b:
                        s_le += calc(f, i)
                table_rows.append(row)
                i += h
                j += 1

            A = calc(f, a) + calc(f, b)
            B = 2 * s_chan
            C = 4 * s_le
            I = (h / 3) * (A + B + C)

            result_simpson = {
                "formatted_f": formatted_f,
                'a':froundf(a),
                'b':froundf(b),
                'n':froundf(n),
                "h": froundf(h),
                "A": froundf(A),
                "B": froundf(B),
                "C": froundf(C),
                "I": froundf(I),
                "h_latex": h_latex,
                "integral_expr": integral_expr,
            }
            active_tab = "Simpson"

        elif 'newton' in request.form:
            # [Logic Newton giữ nguyên]
            try:
                f_input = request.form.get("f_input")
                a = float(request.form.get("a"))
                b = float(request.form.get("b"))
                sll = int(request.form.get("sll"))

                x = symbols('x')
                f = sympify(f_input)
                f_prime = diff(f, x)
                f_prime2 = diff(f_prime, x)

                condition = calc(f, a) * calc(f_prime2, a)
                x0 = a if condition > 0 else b
                condition_value = froundf(condition)

                approximations = []
                fn = x0
                for i in range(sll):
                    fn = fn - calc(f, fn) / calc(f_prime, fn)
                    approximations.append(round(fn, 6))

                if check(approximations, sll):
                    approximations = []
                    fn = x0
                    for i in range(sll):
                        fn = fn - calc(f, fn) / calc(f_prime, fn)
                        approximations.append(round(fn, 9))

                x_n = symbols('x_n')
                f_prime_temp = f_prime.replace(log, lambda arg: Function('ln')(arg))
                newton_formula = latex(x_n - f.subs(x, x_n) / f_prime_temp.subs(x, x_n))

                f_prime_a = abs(calc(f_prime, a))
                f_prime_b = abs(calc(f_prime, b))
                f_prime2_a = abs(calc(f_prime2, a))
                f_prime2_b = abs(calc(f_prime2, b))
                fmin = min(f_prime_a, f_prime_b)
                fmax = max(f_prime2_a, f_prime2_b)
                error_estimation = (fmax / (2 * fmin)) * ((approximations[-1] - approximations[-2]) ** 2)

                sll_minus_1 = sll - 1
                m_over_2m = froundf(fmax / (2 * fmin))

                result_newton = {
                    "formatted_f": f_fomat(f),
                    "formatted_f_prime": f_fomat(f_prime),
                    "formatted_f_prime2": f_fomat(f_prime2),
                    "x0": froundf(x0),
                    "condition": f"{condition_value} {'>' if condition > 0 else '<'} 0",
                    "approximations": approximations,
                    "newton_formula": newton_formula,
                    "f_prime_a": froundf(f_prime_a),
                    "f_prime_b": froundf(f_prime_b),
                    "f_prime2_a": froundf(f_prime2_a),
                    "f_prime2_b": froundf(f_prime2_b),
                    "fmin": froundf(fmin),
                    "fmax": froundf(fmax),
                    "error_estimation": round(error_estimation, 9),
                    "a": froundf(a),
                    "b": froundf(b),
                    "sll": sll,
                    "sll_minus_1": sll_minus_1,
                    "m_over_2m": m_over_2m,
                    "diff_squared": round((approximations[-1] - approximations[-2]) ** 2, 9),
                }
                active_tab = "Newton"
            except Exception as e:
                return render_template("index.html", error=f"Lỗi Newton: {str(e)}", active_tab="Newton")
        elif 'regression' in request.form:
            error = None
            try:
                f_a_input = request.form.get("f_a_input")
                f_b_input = request.form.get("f_b_input")
                f_free_input = request.form.get("f_free_input")
                n = int(request.form.get("n_regression"))
                
                if n <= 0:
                    raise ValueError("Số điểm dữ liệu phải lớn hơn 0.")
                
                arr_x_k = []
                arr_y_k = []
                for i in range(n):
                    x_val = request.form.get(f"x_{i}")
                    y_val = request.form.get(f"y_{i}")
                    if x_val is None or y_val is None or x_val.strip() == "" or y_val.strip() == "":
                        raise ValueError(f"Thiếu giá trị x_{i} hoặc y_{i}")
                    arr_x_k.append(froundf(float(x_val)))
                    arr_y_k.append(froundf(float(y_val)))
                
                x = symbols('x')
                a = symbols('a')
                b = symbols('b')
                y = symbols('y')
                y_k = symbols('y_k')
                
                f_a = sympify(f_a_input)
                f_b = sympify(f_b_input)
                f_free = sympify(f_free_input)
                
                f = a * f_a + b * f_b + f_free
                formatted_f = f_fomat(f)
                
                s_replace = f.subs(x, symbols(f'{x}_k'))
                s = (s_replace - y_k)**2
                s_diff_a = (s_replace - y_k) * f_a.subs(x, symbols(f'{x}_k'))
                s_diff_b = (s_replace - y_k) * f_b.subs(x, symbols(f'{x}_k'))
                
                he_pt1 = [f_a**2, f_a*f_b, f_a*f_free, f_a*y]
                he_pt2 = [f_b*f_a, f_b**2, f_b*f_free, f_b*y]
                
                hpt1 = []
                hpt2 = []

                t_hpt1 = []
                t_hpt2 = []

                for i in range(4):
                    t_hpt1.append(he_pt1[i].subs(x, symbols(f'{x}_k')).subs(y, symbols(f'{y}_k')))
                    t_hpt2.append(he_pt2[i].subs(x, symbols(f'{x}_k')).subs(y, symbols(f'{y}_k')))
                
            
                for i in range(4):
                    hpt1.append(latex(t_hpt1[i].replace(log, lambda arg: Function('ln')(arg))))
                    hpt2.append(latex(t_hpt2[i].replace(log, lambda arg: Function('ln')(arg))))
                

                for i in range(0, 3):
                    temp = tinhtong(he_pt1[i], arr_x_k, n)
                    temp2 = tinhtong(he_pt2[i], arr_x_k, n)
                    he_pt1[i] = round(temp, 6)
                    he_pt2[i] = round(temp2, 6)
                
                t_temp = [0, 0]
                for i in range(0, n):
                    t_temp[0] += calc_two_var(he_pt1[3], arr_x_k[i], arr_y_k[i])
                    t_temp[1] += calc_two_var(he_pt2[3], arr_x_k[i], arr_y_k[i])
                
                he_pt1[3] = round(t_temp[0], 6) * -1
                he_pt2[3] = round(t_temp[1], 6) * -1
                
                f1 = he_pt1[0]*a + he_pt1[1]*b + (he_pt1[2] + he_pt1[3])
                f2 = he_pt2[0]*a + he_pt2[1]*b + (he_pt2[2] + he_pt2[3])
                result = solve_h2(f1, f2)
                
                regression_table = []
                for i in range(n):
                    regression_table.append({
                        'x': arr_x_k[i],
                        'y': arr_y_k[i],
                        'f_a': round(calc(f_a, arr_x_k[i]), 6),
                        'f_b': round(calc(f_b, arr_x_k[i]), 6),
                        'f_free': round(calc(f_free, arr_x_k[i]), 6)
                    })
                

                final_equation = latex(round(result[a],6) * f_a + round(result[b],6) * f_b + f_free).replace('log', 'ln')
                
                
                result_regression = {
                    'formatted_f': formatted_f,
                    'hpt1':hpt1,
                    'hpt2':hpt2,
                    'arr_x_k': arr_x_k,
                    'arr_y_k': arr_y_k,
                    'formatted_s': f_fomat(s),
                    'formatted_s_diff_a': f_fomat(s_diff_a),
                    'formatted_s_diff_b': f_fomat(s_diff_b),
                    'he_pt1': he_pt1,
                    'he_pt2': he_pt2,
                    'f1': f"{he_pt1[0]}a + {he_pt1[1]}b = {-1*(he_pt1[2]+he_pt1[3])}",
                    'f2': f"{he_pt2[0]}a + {he_pt2[1]}b = {-1*(he_pt2[2]+he_pt2[3])}",
                    'a': froundf(result[a]),
                    'b': froundf(result[b]),
                    'n': n,
                    'final_equation': final_equation
                }
                
                active_tab = "Regression"
            except Exception as e:
                error = f"Lỗi Hồi quy: {str(e)}"
                active_tab = "Regression"
        elif 'euler' in request.form:
            # Logic Euler cải tiến
            try:
                y_diff_input = request.form.get("y_diff_input")
                x0 = float(request.form.get("x0"))
                y0 = float(request.form.get("y0"))
                x1 = float(request.form.get("x1"))
                x2 = float(request.form.get("x2"))
                h = float(request.form.get("h"))
                n = int(request.form.get("n"))

                x, y = symbols('x y')
                y_diff = sympify(y_diff_input)
                x0_sym, y0_sym = symbols('x0 y0')
                x1_sym, y1_sym = symbols('x1 y1')
                x2_sym, y2_sym = symbols('x2 y2')

                # Tính y1
                s0_y1 = y0 + h * calc_two_var(y_diff_input, x0, y0)
                a = s0_y1
                sk_y1 = []
                for i in range(n):
                    a = y0 + h / 2 * (calc_two_var(y_diff_input, x0, y0) + calc_two_var(y_diff_input, x1, a))
                    sk_y1.append(round(a, 9))

                y_diff_y1_0 = f_fomat1(y_diff.subs(x, symbols(f'{x}0')).subs(y, symbols(f'{y}0')))
                y_diff_y1_k = f_fomat1(y_diff.subs(x, symbols(f'{x}1')).subs(y, symbols(f'{y}1')))

                # Tính y2
                b = sk_y1[-1]
                s0_y2 = b + h * calc_two_var(y_diff_input, x1, b)
                c = s0_y2
                sk_y2 = []
                for i in range(n):
                    c = b + h / 2 * (calc_two_var(y_diff_input, x1, b) + calc_two_var(y_diff_input, x2, c))
                    sk_y2.append(round(c, 9))

                y_diff_y2_0 = f_fomat1(y_diff.subs(x, symbols(f'{x}1')).subs(y, symbols(f'{y}1')))
                y_diff_y2_k = f_fomat1(y_diff.subs(x, symbols(f'{x}2')).subs(y, symbols(f'{y}2')))

                result_euler = {
                    "y_diff": f_fomat1(y_diff_input),
                    "x0": froundf(x0),
                    "y0": froundf(y0),
                    "x1": froundf(x1),
                    "x2": froundf(x2),
                    "h": froundf(h),
                    "n": n,
                    "s0_y1": froundf(s0_y1),
                    "y_diff_y1_0": y_diff_y1_0,
                    "y_diff_y1_k": y_diff_y1_k,
                    "sk_y1": sk_y1,
                    "s0_y2": froundf(s0_y2),
                    "y_diff_y2_0": y_diff_y2_0,
                    "y_diff_y2_k": y_diff_y2_k,
                    "sk_y2": sk_y2,
                }
                active_tab = "Euler"

            except ValueError:
                error = "Lỗi Euler: Vui lòng nhập các giá trị số hợp lệ."
                active_tab = "Euler"
            except SympifyError:
                error = "Lỗi Euler: Biểu thức y' không hợp lệ."
                active_tab = "Euler"
            except Exception as e:
                error = f"Lỗi Euler: {str(e)}"
                active_tab = "Euler"

    return render_template(
        "index.html",
        result_simpson=result_simpson,
        result_newton=result_newton,
        result_regression=result_regression, 
        result_euler=result_euler,
        table_rows=table_rows,
        regression_table=regression_table, 
        active_tab=active_tab,
        error=error  # Thêm error vào template
    )


if __name__ == "__main__":
    app.run(debug=True)