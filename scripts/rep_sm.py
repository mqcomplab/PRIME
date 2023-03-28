import sys
sys.path.insert(0, '../')
import modules as mod

mod.gen_method_max(trim=None, n_ary="SM")
mod.gen_method_max(trim=0.1, n_ary="SM")
mod.gen_method_max(trim=0.2, n_ary="SM")
mod.gen_method_max(trim=0.3, n_ary="SM")

