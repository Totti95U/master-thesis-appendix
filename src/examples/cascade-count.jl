
marker = ("*",) # empty shift
# marker = (2, "*", 2, 1, 1) # γ_q (same as Lipa's B*BAA)
# marker = (2, 1, "*", 2, 1, 1)
# marker = (2, 1, "*", 2, 1) # γ_r
# marker = (1, "*", 2, 1) # γ_s
# marker = [(1, "*", 2, 1), (2, 2, 1, "*", 2, 2, 2, 1)]
# marker = [(2, 1, "*", 2, 1), (1, 1, 1, "*", 2, 1), (2, 1, "*", 2, 2, 2)]
# marker = [(1, "*", 2, 1), (2, 1, "*", 2, 2, 2)]
# marker = (1, 1, 2, "*", 2, 1, 1) # γ_p
# marker = () # full shift

# marker = (2, 2, "*", 2, 1, 1) # Lipa's BB*BAA
# marker = (1, 2, "*", 2, 1, 1) # Lipa's AB*BAA
# markers = [(1, "*", 2, 1, 1), (1, "*", 2, 1, 2, 2, 1)] # Lipa's A*BAA and A*BABBA

marker1 = (2, 1, "*", 2, 1, 1) # こっちが superset
marker2 = (2, 1, "*", 2, 1) # こっちが subset

res = gcd2_cascade(marker2, 16) .- gcd2_cascade(marker1, 16)

println("-" ^ 60)
# println("count_cascade:")
# display(count_cascade(marker, 16)')
println("gcd2_cascade: $marker1 \\ $marker2")
display(res')
# println("even1_orbit:")
# display(even1_cascade(marker, 16)')

"""
Remark: Y1\\Y2 の 2-cascade 分解について計算したい場合は Y2 に関する 2-cascade 分解の結果から Y1 に関する 2-cascade 分解の結果を引けば良い
※ X > Y > Z のとき (X\\Y)\\(X\\Z) = Y\\Z であることを利用
""";
