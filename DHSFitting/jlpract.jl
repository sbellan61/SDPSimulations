using Distributions
using Devectorize

K = 100
s__L = ones(K)

p_m_bef = zeros(K)
p_f_bef = zeros(K)


p_m_bef = rand(K)
p_f_bef = rand(K)

states = rand(K,9)
states = DataFrame(0, K, 9)

states = DataFrame(0, K, 9)
colnames!(states, ["s__", "mb_a1", "mb_a2", "mb_", "f_ba1", "f_ba2", "f_b", "hb1b2", "hb2b1"])
states["s__"] = ones(K)
statesL = states

function pre_couple(states, p_m_bef, p_f_bef)
    s__   = s__L .* (1-p_m_bef) .* (1-p_f_bef)
    mb_a1 = s__L .* p_m_bef .* (1-p_f_bef)
    mb_a2 = mb_a1L .* (1 - p_f_bef)
    mb_   = mb_a2L .* (1 - p_f_bef) + mb_L .* (1 - p_f_bef)
    f_ba1 = s__L .* p_f_bef .* (1-p_m_bef)
    f_ba2 = f_ba1L .* (1 - p_m_bef)
    f_b   = f_ba2L .* (1 - p_m_bef) + f_bL .* (1 - p_m_bef)
    p_mfirst = p_m_bef / (p_m_bef+p_f_bef)
    p_ffirst = 1-p_mfirst
    p_mfirst[is_na(p_mfirst)] = 0
    p_ffirst[is_na(p_ffirst)] = 0                
    hb1b2 = hb1b2L + p_mfirst  .*  s__L .* p_m_bef .* p_f_bef + (mb_a1L + mb_a2L + mb_L)  .*  p_f_bef
    hb2b1 = hb2b1L + p_ffirst  .*  s__L .* p_m_bef .* p_f_bef + (f_ba1L + f_ba2L + f_bL)  .*  p_m_bef
end

                    

function dynupdate(states, logic)
    with(states[logic,:], {
                           x1 = x1 .* 2
                           x2 = par1 .* x3 .*x 2
                           x3 = par3 .* (x1 + x3)^2})
    return(states)
end



   

    
