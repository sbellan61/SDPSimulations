(s__, mb_a1, mb_a2, mb_, f_ba1, f_ba2, f_b, hb1b2, hb2b1) = 1:9
function pre_coupleMat(serostates, sexually_active) 
    temp = serostates[sexually_active,:]
    serostates[sexually_active,s__]   = temp[:,s__] .* (1-p_m_bef) .* (1-p_f_bef)
    serostates[sexually_active,mb_a1] = temp[:,s__] .* p_m_bef .* (1-p_f_bef)
    serostates[sexually_active,mb_a2] = temp[:,mb_a1] .* (1 - p_f_bef)
    serostates[sexually_active,mb_] = temp[:,mb_a2] .* (1 - p_f_bef) + temp[:,mb_] .* (1 - p_f_bef)
    serostates[sexually_active,f_ba1] = temp[:,s__] .* p_f_bef .* (1-p_m_bef)
    serostates[sexually_active,f_ba2] = temp[:,f_ba1] .* (1 - p_m_bef)
    serostates[sexually_active,f_b] = temp[:,f_ba2] .* (1 - p_m_bef) + temp[:,f_b] .* (1 - p_m_bef)
    serostates[sexually_active,hb1b2] = temp[:,hb1b2] + .5  .*  temp[:,s__] .* p_m_bef .* p_f_bef + (temp[:,mb_a1] + temp[:,mb_a2] + temp[:,mb_])  .*  p_f_bef
    serostates[sexually_active,hb2b1] = temp[:,hb2b1] + .5  .*  temp[:,s__] .* p_m_bef .* p_f_bef + (temp[:,f_ba1] + temp[:,f_ba2] + temp[:,f_b])  .*  p_m_bef
    return serostates
end

function it_serostates(serostates, its)
    for ii=[1:its]
        serostates = pre_coupleMat(serostates, sexually_activeMat[:,1])
    end
    return(serostates)
end

n = 10^4
k = 9
serostates = zeros(n,k)
serostates[:,s__] = 1
serostates
sexually_activeMat = randbool(n, its)
sexually_activeMat[:,1]
serostates[sexually_activeMat[:,1],s__]

p_m_bef = .012
p_f_bef = .07
@time it_serostates(serostates, 12*40)
