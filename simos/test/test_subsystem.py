# import simos
# import pytest
# import qutip as qu
# import sympy as sp
# import numpy as np

# @pytest.mark.parametrize("method", ['qutip','numpy','sympy','sparse' ])
# class TestSubsystem:
#     def test_subsystem_1(self,method):
#         A = {'name':'A', 'val':1}
#         B = {'name':'B', 'val':1/2}
#         C = {'name':'C','val': 3/2}
#         s = simos.create_system([A, B, C] ,method=method)
#         test,left, reverse = simos.subsystem(s, s.id, "B", keep  = True)
#         test2 = simos.reverse_subsystem(s, test, left, reverse)

#         def _all(test):
#             if method == 'qutip':
#                 return test
#             elif method == 'numpy' or method == 'numba':
#                 return test.all()
#             elif method == 'sympy':
#                 return test
#             elif method == 'sparse':
#                 return test.toarray().all()

#         assert _all((s.id == test2))

#         return None 

#     def test_subsystem_2(self,method):
#         A = {'name':'A', 'val':1}
#         B = {'name':'B', 'val':1/2}
#         C = {'name':'C','val': 3/2}
#         s = simos.create_system((A, B, C) ,method=method)
#         test,left, reverse = simos.subsystem(s, s.id, "B", keep  = True)
#         test2 = simos.reverse_subsystem(s, test, left, reverse)

#         def _all(test):
#             if method == 'qutip':
#                 return test
#             elif method == 'numpy' or method == 'numba':
#                 return test.all()
#             elif method == 'sympy':
#                 return test
#             elif method == 'sparse':
#                 return test.toarray().all()

#         assert _all((s.id == test2))

#         return None         

#     def test_subsystem_3(self,method):
#         A = {'val': 1/2, 'name':'A'}  
#         B = {'val': 1, 'name':'B'} 
#         C = {'val': 3/2, 'name':'C'} 
#         GS= {'val': 0, 'name':'GS'}  
#         ES = {'val':0, 'name':'ES'}
#         s = simos.System([([GS, B, C],ES), A], method = method)
#         test,left, reverse = simos.subsystem(s, s.id, "B", keep  = True)
#         test2 = simos.reverse_subsystem(s, test, left, reverse)

#         def _all(test):
#             if method == 'qutip':
#                 return test
#             elif method == 'numpy' or method == 'numba':
#                 return test.all()
#             elif method == 'sympy':
#                 return test
#             elif method == 'sparse':
#                 return test.toarray().all()

#         assert _all((s.id == test2))

#         return None     
    

#     def test_subsystem_4(self,method):
#         A = {'val': 1/2, 'name':'A'}  
#         B = {'val': 1, 'name':'B'} 
#         C = {'val': 3/2, 'name':'C'} 
#         GS= {'val': 0, 'name':'GS'}  
#         ES = {'val':0, 'name':'ES'}
#         s = simos.System([([GS, B, C],ES), A], method = method)
#         test,left, reverse = simos.subsystem(s, s.id, ["GS", "ES"], keep  = True)
#         test2 = simos.reverse_subsystem(s, test, left, reverse)

#         def _all(test):
#             if method == 'qutip':
#                 return test
#             elif method == 'numpy' or method == 'numba':
#                 return test.all()
#             elif method == 'sympy':
#                 return test
#             elif method == 'sparse':
#                 return test.toarray().all()

#         assert _all((s.id == test2))

#         return None  
    
    # this should throw an error verify that
    # def test_subsystem_5(self,method):
    #     A = {'val': 1/2, 'name':'A'}  
    #     B = {'val': 1, 'name':'B'} 
    #     C = {'val': 3/2, 'name':'C'} 
    #     GS= {'val': 0, 'name':'GS'}  
    #     ES = {'val':0, 'name':'ES'}
    #     s = simos.System([([GS, B, C],ES), A], method = method)
    #     test,left, reverse = simos.subsystem(s, s.id, ["GS", "ES", "A"], keep  = True)
    #     test2 = simos.reverse_subsystem(s, test, left, reverse)

    #     def _all(test):
    #         if method == 'qutip':
    #             return test
    #         elif method == 'numpy' or method == 'numba':
    #             return test.all()
    #         elif method == 'sympy':
    #             return test
    #         elif method == 'sparse':
    #             return test.toarray().all()

    #     assert _all((s.id == test2))

    #     return None  
    
    # this should throw an error verify that
    # def test_subsystem_6(self,method):
    #     A = {'val': 1/2, 'name':'A'}  
    #     B = {'val': 1, 'name':'B'} 
    #     C = {'val': 3/2, 'name':'C'} 
    #     GS= {'val': 0, 'name':'GS'}  
    #     ES = {'val':0, 'name':'ES'}
    #     s = simos.System([([GS, B, C],ES), A], method = method)
    #     test,left, reverse = simos.subsystem(s, s.id, ["GS", "C", "A"], keep  = True)
    #     test2 = simos.reverse_subsystem(s, test, left, reverse)

    #     def _all(test):
    #         if method == 'qutip':
    #             return test
    #         elif method == 'numpy' or method == 'numba':
    #             return test.all()
    #         elif method == 'sympy':
    #             return test
    #         elif method == 'sparse':
    #             return test.toarray().all()

    #     assert _all((s.id == test2))

    #     return None  