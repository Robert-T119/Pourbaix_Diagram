from django.urls import path
from . import views


urlpatterns = [
    path('', views.nickelca, name='nickelca')
]
